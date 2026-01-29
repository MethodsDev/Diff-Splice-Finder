# Control Groups Feature - Implementation Summary

## Overview
Added functionality to specify control sample groups that all other sample types are compared against, rather than performing all pairwise comparisons.

## Changes Made

### 1. Main Pipeline Script (`run_diff_splice_analysis.py`)

#### New Parameter
- `--control_groups`: Comma-separated list of control group names
  - Example: `--control_groups control`
  - Example: `--control_groups control,wildtype`
  - Mutually exclusive with `--contrast` parameter

#### Logic Changes
- Modified contrast generation in `run_edgeR()` function (lines ~350-385):
  - If `--control_groups` is specified:
    - Parses control groups from comma-separated string
    - Validates that all specified control groups exist in the data
    - Generates contrasts comparing each non-control group vs control(s)
    - Format: "TreatmentGroup-ControlGroup1,ControlGroup2"
  - If not specified:
    - Uses original behavior (all pairwise comparisons or single contrast)

#### Validation
- Added validation to prevent simultaneous use of `--contrast` and `--control_groups`
- Clear error messages if control groups are not found in the data

### 2. R Analysis Script (`util/run_edgeR_analysis.R`)

#### Updated Contrast Parsing (lines 135-194)
- Modified to handle comma-separated control groups in contrast strings
- Examples:
  - "TDP43-control" → single control
  - "TDP43-control,wildtype" → multiple controls pooled together
- For multiple controls:
  - Weights are distributed equally: each control gets `-1/n_controls`
  - This averages the control groups in the statistical comparison
- Improved error messages for validation

### 3. Documentation

#### Examples Created
- `examples/run_with_control_groups.sh`: Shell script demonstrating usage
- `examples/sample_metadata_with_controls.tsv`: Example metadata with multiple groups

#### Updated Files
- `README.md`: Added "With Control Groups" section showing examples
- `examples/PARAMETER_GUIDE.md`: 
  - Added detailed description of `--control_groups` parameter
  - Added section on control-based comparison workflows
  - Documented benefits (reduced multiple testing, more intuitive results)

### 4. Testing
- `testing/test_metadata_control.tsv`: Test metadata file created
- Validated functionality with real data:
  - Successfully identified control and treatment groups
  - Correctly generated contrasts
  - Produced expected output files
  - Error handling works correctly

## Usage Examples

### Single Control Group
```bash
./run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata.tsv \
    --output_dir output_with_controls \
    --control_groups control \
    --min_delta_psi 0.1 \
    --fdr_threshold 0.05
```

With groups: TDP43, FUS, TARDBP, control
- This generates 3 comparisons:
  - TDP43 vs control
  - FUS vs control
  - TARDBP vs control

### Multiple Control Groups
```bash
./run_diff_splice_analysis.py \
    --matrix data/intron_counts.matrix \
    --samples examples/sample_metadata_with_controls.tsv \
    --output_dir output_pooled_controls \
    --control_groups control,wildtype
```

With groups: TDP43, FUS, TARDBP, control, wildtype
- This generates 3 comparisons:
  - TDP43 vs control,wildtype (pooled)
  - FUS vs control,wildtype (pooled)
  - TARDBP vs control,wildtype (pooled)

## Benefits

1. **Focused Analysis**: Compare only treatment vs control (biologically relevant)
2. **Reduced Multiple Testing**: Fewer comparisons → better FDR control
3. **Statistical Power**: Pooling multiple control types increases sample size
4. **Interpretability**: Results directly answer "what changed vs control?"
5. **Flexibility**: Can specify one or many control groups

## Technical Details

### Statistical Approach
When multiple control groups are specified, they are pooled by:
- Creating a contrast vector where each control group gets weight `-1/n_controls`
- This effectively tests: Treatment vs Average(Controls)
- Maintains proper statistical properties in the GLM framework

### Contrast Format in Output
- Single control: `TDP43_vs_control`
- Multiple controls: `TDP43_vs_control_wildtype`
- This appears in the `contrast` column of results files

## Backward Compatibility

All existing functionality is preserved:
- Default behavior unchanged (all pairwise comparisons)
- `--contrast` parameter still works as before
- No changes to output file formats or column names
- Existing scripts and workflows continue to work

## Files Modified

1. `run_diff_splice_analysis.py`
2. `util/run_edgeR_analysis.R`
3. `README.md`
4. `examples/PARAMETER_GUIDE.md`
5. `examples/run_with_control_groups.sh` (new)
6. `examples/sample_metadata_with_controls.tsv` (new)
7. `testing/test_metadata_control.tsv` (new)
