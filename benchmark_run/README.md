# Benchmark Usage Guide

## Prerequisites
Ensure that benchmark source files are generated in the current directory before proceeding.

## Running Benchmarks

### Step 1: Configure the Benchmark Runner
Edit `json_runner.sh` to specify:
- The binary source to benchmark
- The output directory for results

### Step 2: Execute the Benchmark
```bash
bash json_runner.sh
```

### Step 3: Results
After execution, a JSON file containing benchmark results will be generated in the specified output directory.

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

### Step 2: Configure Data Processing
Open `processdata.py` and modify the main function to include:
- Path to your benchmark JSON data file
- Desired output directory for generated plots

### Step 3: Generate Plots
```bash
python3 processdata.py
```

### Step 4: Results
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

## File Structure
```
.
├── json_runner.sh          # Benchmark execution script
├── processdata.py          #Data processing and visualization
├── requirements.txt        # Python dependencies
└── myenv/                  # Python virtual environment (created)
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
