# Benchmark Usage Guide

## Prerequisites
Ensure that data files are generated in the server directory before proceeding.

## Data Processing and Visualization

### Step 1: Set Up Python Environment
```bash
# Create a virtual environment
python3 -m venv myenv

# Activate the virtual environment
source myenv/bin/activate

# Install required dependencies
pip3 install -r requirements.txt
```

### Step 2: Generate Plots
```bash
run_python.sh
```

### Step 3: Results
Visualization plots will be generated in the specified output directory.

## Environment Management

### Activating the Virtual Environment
Whenever you need to work with the benchmarking tools:
```bash
source myenv/bin/activate
```

### Deactivating the Virtual Environment
When finished:
```bash
deactivate
```

## Troubleshooting

### Virtual Environment Issues
If you encounter permission errors:
```bash
# Recreate the virtual environment if needed
rm -rf myenv
python3 -m venv myenv
```

### Dependency Installation Issues
If pip fails to install requirements:
```bash
# Upgrade pip first
pip3 install --upgrade pip

# Then install requirements
pip3 install -r requirements.txt
```

### Benchmark Execution Issues
Ensure that:
- Benchmark binaries are properly compiled
- Output directory has write permissions
- JSON runner script has execute permissions: `chmod +x json_runner.sh`
