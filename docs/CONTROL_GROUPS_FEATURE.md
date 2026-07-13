# Control Groups Feature

## Historical Note

This document describes a previous multi-contrast/control-group orchestration
feature. The current `DSF.py` interface requires one
explicit contrast per invocation:

```bash
--contrast GroupA,GroupB
```

Pooled semicolon-delimited controls and `--control_groups` are not part of the
current main pipeline contract. To compare several treatment groups against a
control, run the pipeline once per treatment/control contrast and use separate
output directories.

The active contrast behavior is documented in `README.md` and
`examples/PARAMETER_GUIDE.md`.
