# Plasma Downstream Analysis Modules

## Status: Coming Soon

This directory will contain plasma-specific downstream analysis modules.

## Expected Structure

When plasma downstream scripts are provided, they will be integrated here following the same pattern as saliva modules:

```
downstream_plasma/
├── analysis_1.nf
├── analysis_2.nf
├── ...
└── aggregate.nf
```

## Adding Plasma Modules

1. Create `.nf` module files in this directory
2. Update `main.nf` to include the new modules
3. Add plasma-specific downstream flags (e.g., `--run_plasma_analysis_1`)
4. Update the `RUN_DOWNSTREAM` workflow to handle plasma mode

## Contact

Provide plasma downstream scripts to extend this pipeline.
