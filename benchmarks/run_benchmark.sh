#!/bin/bash
# CarrierCapture.jl vs CarrierCapture.py Benchmark Runner
# =======================================================
#
# This script runs the full benchmark to validate that CarrierCapture.py
# gives identical results to CarrierCapture.jl for the Sn_Zn in ZnO example.
#
# Requirements:
#   - Julia with CarrierCapture.jl installed
#   - Python with CarrierCapture.py installed
#
# Usage:
#   ./benchmarks/run_benchmark.sh

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "============================================================"
echo "CarrierCapture.jl vs CarrierCapture.py Benchmark"
echo "============================================================"
echo ""

# Check if we're in the right directory
if [ ! -f "benchmarks/run_julia_reference.jl" ]; then
    echo -e "${RED}✗ ERROR${NC}: Must run from repository root directory"
    echo "Usage: ./benchmarks/run_benchmark.sh"
    exit 1
fi

# Create directories
echo "Creating output directories..."
mkdir -p benchmarks/reference_data
mkdir -p benchmarks/results
echo -e "${GREEN}✓${NC} Directories ready"
echo ""

# Step 1: Run Julia reference
echo "============================================================"
echo "Step 1: Running Julia reference calculation..."
echo "============================================================"
echo ""

if ! command -v julia &> /dev/null; then
    echo -e "${RED}✗ ERROR${NC}: Julia not found in PATH"
    echo "Please install Julia or add it to your PATH"
    exit 1
fi

echo "Julia version:"
julia --version
echo ""

echo "Running Julia script..."
if julia benchmarks/run_julia_reference.jl; then
    echo ""
    echo -e "${GREEN}✓ Julia reference complete${NC}"
else
    echo ""
    echo -e "${RED}✗ Julia reference failed${NC}"
    echo "Check error messages above"
    exit 1
fi

echo ""

# Step 2: Run Python benchmark
echo "============================================================"
echo "Step 2: Running Python benchmark and comparison..."
echo "============================================================"
echo ""

if ! command -v python &> /dev/null; then
    echo -e "${RED}✗ ERROR${NC}: Python not found in PATH"
    echo "Please install Python or add it to your PATH"
    exit 1
fi

echo "Python version:"
python --version
echo ""

echo "Running Python benchmark..."
if python benchmarks/benchmark_sn_zn.py; then
    echo ""
    echo -e "${GREEN}✓ Benchmark complete - ALL TESTS PASSED${NC}"
    BENCHMARK_RESULT=0
else
    echo ""
    echo -e "${YELLOW}⚠ Benchmark complete - SOME TESTS FAILED${NC}"
    echo "Check benchmark report for details"
    BENCHMARK_RESULT=1
fi

echo ""

# Summary
echo "============================================================"
echo "Benchmark Summary"
echo "============================================================"
echo ""
echo "Results saved to:"
echo "  - Julia reference:  benchmarks/reference_data/sn_zn_julia_reference.json"
echo "  - Benchmark report: benchmarks/results/benchmark_report.json"
echo ""

if [ $BENCHMARK_RESULT -eq 0 ]; then
    echo -e "${GREEN}✓ Success!${NC} Python implementation validated against Julia reference."
    echo ""
    echo "You can now add these results to the README:"
    echo "  - Update the 'Performance' section with actual benchmark data"
    echo "  - Add a 'Validation' section showing test pass/fail status"
else
    echo -e "${YELLOW}⚠ Warning${NC}: Some validation tests failed."
    echo "Review benchmarks/results/benchmark_report.json for details."
fi

echo ""
echo "============================================================"

exit $BENCHMARK_RESULT
