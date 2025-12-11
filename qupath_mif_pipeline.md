# QuPath mIF pipeline for complex TME (Tregs and T cells)

This document outlines a practical QuPath automation plan for multiplex immunofluorescence (mIF) slides with T cell subsets and Tregs (e.g., CD3, CD8, CD4, FOXP3, PD-1/PD-L1, pan-cytokeratin). It includes configuration guidance, quality controls, and a Groovy batch script template you can adapt for your project.

## 1) Channel and project setup
- Standardize channel names in the project (e.g., `DAPI`, `CD3`, `CD8`, `CD4`, `FOXP3`, `PD1`, `PDL1`, `PanCK`).
- Ensure pixel calibration is correct for the scanner; confirm the unit is microns.
- Add slides to a QuPath project and store versioned parameter files (YAML/JSON) that capture thresholds and classifier names.
- (Optional) Import ROI annotations (GeoJSON/JSON) for tumor regions or pathologist-marked hot/cold spots.

## 2) Pre-processing
- Perform background subtraction per channel (rolling ball or flat-field) and correct illumination if a shading profile exists.
- For crosstalk-heavy panels, apply spectral unmixing results before running cell detection.
- Normalize intensities per channel to robust percentiles (e.g., 1st–99th) to stabilize downstream thresholds.

## 3) Tissue and compartment masks
- Train/use a pixel classifier to segment **tissue vs. glass** and compartments (e.g., tumor vs. stroma) using PanCK and morphology cues.
- Export masks as probability maps or QuPath objects; keep classifier names versioned.
- Downstream thresholds should be compartment-aware (e.g., PanCK-positive tumor nests vs. stroma).

## 4) Cell detection (nuclear-driven)
- Use **Watershed Cell Detection** on DAPI with typical mIF nuclear radii (e.g., nucleus radius 5–8 µm, background radius 8–12 µm).
- Enable splitting sensitivity to avoid clusters in lymphocyte-dense regions.
- Expand cytoplasm measurement (e.g., 2–3 µm) for membrane markers like PD-1/PD-L1.

## 5) Marker gating and phenotypes
- Derive thresholds from control slides or percentile-based rules; lock them in a YAML/JSON file (example below).
- Define base marker positivity per compartment (tumor/stroma), then composite phenotypes:
  - **CD3+CD8+** (cytotoxic T cells)
  - **CD3+CD4+FOXP3+** (Tregs)
  - **CD3+CD4+FOXP3-** (helper T cells)
  - **CD3+PD1+** (exhausted/activated T cells)
  - **PanCK+PDL1+** (PD-L1+ tumor)
  - **CD68+PDL1+** (myeloid PD-L1; include if macrophage channel is present)
- Store thresholds like:
  ```yaml
  thresholds:
    CD3: 1200
    CD8: 950
    CD4: 900
    FOXP3: 800
    PD1: 700
    PDL1: 750
    PanCK: 850
  classifiers:
    tissue: "Tissue_vs_Glass_v2"
    compartments: "Tumor_Stroma_v1"
  ```

## 6) Quality control
- Export per-slide QC grids (PNG) with overlays for tissue mask, detections, and phenotypes.
- Log counts/densities per compartment and % positivity per marker; flag outliers (e.g., >3 MAD from cohort median).
- Keep a manifest tying slide IDs to parameter file versions and ROI names.

## 7) Batch export
- Export **per-cell tables** (CSV/TSV) with: slide ID, ROI, x/y in microns, nucleus/cytoplasm measurements, marker intensities, phenotype labels, and compartment labels.
- Optionally export GeoJSON of detections and probability maps for spatial analysis in Python/R.

## 8) Groovy script template (run via `qupath script`)
Below is a concise Groovy template to batch process all images in a project. Update classifier names, thresholds, and output paths before running.

```groovy
// QuPath Groovy batch for mIF TME with Tregs/T cells
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.objects.PathObject

// ---- User-configurable parameters ----
// Thresholds should match your YAML/JSON records
thresholds = [
  CD3: 1200,
  CD8: 950,
  CD4: 900,
  FOXP3: 800,
  PD1: 700,
  PDL1: 750,
  PanCK: 850
]

def tissueClassifierName = "Tissue_vs_Glass_v2"
def compartmentClassifierName = "Tumor_Stroma_v1"
def exportDir = buildFilePath(PROJECT_BASE_DIR, "exports", "cells")
mkdirs(exportDir)

// ---- Helper: apply phenotype from thresholds ----
PathClass cd3 = PathClassFactory.getPathClass("CD3+")
PathClass cd8 = PathClassFactory.getPathClass("CD8+")
PathClass cd4 = PathClassFactory.getPathClass("CD4+")
PathClass foxp3 = PathClassFactory.getPathClass("FOXP3+")
PathClass pd1 = PathClassFactory.getPathClass("PD1+")
PathClass pdl1 = PathClassFactory.getPathClass("PDL1+")
PathClass panck = PathClassFactory.getPathClass("PanCK+")

// Composite classes
PathClass treg = PathClassFactory.getPathClass("Treg (CD3+CD4+FOXP3+)")
PathClass cd8eff = PathClassFactory.getPathClass("CD8 T cell (CD3+CD8+)")
PathClass helper = PathClassFactory.getPathClass("CD4 helper (CD3+CD4+FOXP3-)")
PathClass cd3pd1 = PathClassFactory.getPathClass("CD3+PD1+")
PathClass tumorPDL1 = PathClassFactory.getPathClass("Tumor PDL1+ (PanCK+)")

// ---- Processing loop ----
for (imageData in getProject().getImageList()) {
  print "Processing ${imageData.getImageName()}"
  setBatchProjectAndImage(imageData)

  // Run tissue & compartment classifiers
  runPixelClassifier(tissueClassifierName)
  runPixelClassifier(compartmentClassifierName)

  // Cell detection tuned for lymphocyte-rich TME
  selectAnnotations()
  runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', [
    'detectionImage'                : 'DAPI',
    'requestedPixelSizeMicrons'     : 0.5,
    'backgroundRadiusMicrons'       : 10.0,
    'medianRadiusMicrons'           : 0.0,
    'sigmaMicrons'                  : 1.5,
    'minAreaMicrons'                : 20.0,
    'maxAreaMicrons'                : 400.0,
    'threshold'                     : 80.0,
    'watershedPostProcess'          : true,
    'cellExpansionMicrons'          : 2.5,
    'includeNuclei'                 : true,
    'smoothBoundaries'              : true,
    'makeMeasurements'              : true
  ])

  // Assign base marker positivity
  def cells = getCellObjects()
  cells.each { cell ->
    def m = cell.getMeasurementList()
    def classes = [] as Set
    if (m.getMeasurementValue('Mean CD3')   > thresholds.CD3)   classes << cd3
    if (m.getMeasurementValue('Mean CD8')   > thresholds.CD8)   classes << cd8
    if (m.getMeasurementValue('Mean CD4')   > thresholds.CD4)   classes << cd4
    if (m.getMeasurementValue('Mean FOXP3') > thresholds.FOXP3) classes << foxp3
    if (m.getMeasurementValue('Mean PD1')   > thresholds.PD1)   classes << pd1
    if (m.getMeasurementValue('Mean PDL1')  > thresholds.PDL1)  classes << pdl1
    if (m.getMeasurementValue('Mean PanCK') > thresholds.PanCK) classes << panck

    // Composite phenotypes
    if (classes.contains(cd3) && classes.contains(cd8)) cell.setPathClass(cd8eff)
    else if (classes.contains(cd3) && classes.contains(cd4) && classes.contains(foxp3)) cell.setPathClass(treg)
    else if (classes.contains(cd3) && classes.contains(cd4) && !classes.contains(foxp3)) cell.setPathClass(helper)
    else if (classes.contains(cd3) && classes.contains(pd1)) cell.setPathClass(cd3pd1)
    else if (classes.contains(panck) && classes.contains(pdl1)) cell.setPathClass(tumorPDL1)
    else if (!classes.isEmpty()) cell.setPathClass(classes.iterator().next()) // keep single marker label
    else cell.setPathClass(null)
  }
  fireHierarchyUpdate()

  // Export per-cell table
  def out = buildFilePath(exportDir, imageData.getImageName() + '.csv')
  exportCellMeasurements(out)
}
println 'Batch complete.'
```

### Notes on adapting the script
- Replace `Mean <marker>` with your measurement names (check the exact feature labels in QuPath after running detection).
- If you gate on **compartments**, add checks for the assigned parent annotation class (tumor vs. stroma) before applying composite phenotypes.
- Tune `threshold` and area settings in the detection block for your scanner resolution.
- Use `runClassifier` outputs to filter ROIs: e.g., `selectObjectsByClassification("Tumor")` before detection if you only want tumor regions.

## 9) Downstream analysis hooks
- After exporting CSVs, use Python/R to compute densities per compartment, co-expression frequencies, and spatial proximity (e.g., kNN within 50 µm for CD8→PanCK or Treg→CD8 interactions).
- Keep all parameter files, classifier versions, and export manifests under version control for auditability.

## 10) Reproducibility checklist
- Freeze software versions (QuPath, Java runtime, classifier model versions).
- Track threshold YAML/JSON alongside export outputs.
- Keep a minimal test slide to validate the pipeline when upgrading QuPath or classifiers.
