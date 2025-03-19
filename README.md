# Funnel and Pnear Plotter for Helical Bundle Predictions and FastRelax

## Description
This script generates funnel plots and calculates Pnear values from `helical_bundle_predict` and `fastrelax` near-native log files. It is designed to visualize the energy landscape and assess near-native structure quality.

## Usage
```bash
python plot_FR+HBP_funnel.py <any_string>
```

### Requirements
- Run the script in the directory containing `bundleout.log`.
- The `score.sc` file must be present in the `FRfill/` subdirectory.

## Directory Structure
```
RUN_DIRECTORY/
├── bundleout.log
└── FRfill/
    └── score.sc
```

## Output
The script produces:
- A funnel plot visualizing energy vs. RMSD.
- Pnear calculations to quantify near-native convergence.

## Dependencies
- Python 3+
- Matplotlib
- Seaborn
- Scipy
- NumPy
- Pandas

## Author
Allon Goldberg  
Research Assistant  
Flatiron Institute  
2/2025

