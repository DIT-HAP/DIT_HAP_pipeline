# Guide
!!! CRITICAL INSTRUCTIONS FOR THE AI AGENT !!!
1.  **STYLE AND STRUCTURE TEMPLATE ONLY**: Use this template exclusively for its:
    *   Import organization
    *   Logging setup pattern
    *   Use of Pydantic models for validation
    *   Function decorators (like `@logger.catch`)
    *   Overall code structure and formatting
2.  **IGNORE AND REPLACE EXAMPLE CONTENT**: You MUST replace all example content, including:
    *   The specific function names (`example_placeholder_function`, `another_placeholder_function`, `set_args`) and their logic.
    *   The specific docstring content.
3.  **CREATE ORIGINAL CODE**: Your task is to create **ORIGINAL CODE** for the user's specific request, but to format and structure it according to the patterns you see here.
4.  **DOCSTRING STYLE GUIDELINES**:
    *   **IF** the user explicitly requests documentation for a function/method/class that is **public-facing or intended for external use**, use **NumPy style** docstrings.
    *   **OTHERWISE**, use a simple, concise docstring style similar to the existing template (e.g., one-liner or brief description without sections like Parameters, Returns, Examples).
!!! END CRITICAL INSTRUCTIONS !!!

-------------------------------------------------------------------------------

# Template
```python
"""
The python script template for the DIT-HAP project.

This is just a template and you do not need to keep the same document content. You can replace the `[Content in the brackets]` with the proper content.

Typical Usage:
    python [script_name].py --input [input_file] --output [output_file]

Input: [Input file]
Output: [Output file]
[Other information]
"""

# =============================== Imports ===============================
# import necessary libraries, just keep the ones you need, the following are just examples
import sys
import argparse
from pathlib import Path
from loguru import logger
from typing import List, Optional, Dict, Tuple
from pydantic import BaseModel, Field, field_validator
import numpy as np
import pandas as pd

# The following is for plotting, you can remove it if not needed
# from matplotlib import pyplot as plt


# =============================== Constants ===============================
# The following is for plotting, you can remove it if not needed
# plt.style.use('path/to/style.mplstyle')
# AX_WIDTH, AX_HEIGHT = plt.rcParams['figure.figsize']
# COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']


# =============================== Configuration & Models ===============================
# REPLACE THIS ENTIRE CLASS WITH MODELS RELEVANT TO YOUR TASK.
# THIS IS JUST AN EXAMPLE OF HOW TO STRUCTURE A VALIDATION MODEL.
class InputOutputConfig(BaseModel):
    """Pydantic model for validating and managing input/output paths."""
    input_file: Path = Field(..., description="Path to the input file")
    output_file: Path = Field(..., description="Path to the output file")

    @field_validator('input_file')
    def validate_input_file(cls, v):
        if not v.exists():
            raise ValueError(f"Input file does not exist: {v}")
        return v
    
    @field_validator('output_file')
    def validate_output_file(cls, v):
        v.parent.mkdir(parents=True, exist_ok=True) # Create dir if it doesn't exist
        return v
    
    class Config:
        frozen = True # Makes the model immutable after creation

# REPLACE THIS ENTIRE CLASS WITH A MODEL FOR YOUR SPECIFIC RESULTS.
class AnalysisResult(BaseModel):
    """Pydantic model to hold and validate the results of the analysis."""
    total_items_processed: int = Field(..., ge=0, description="Total number of items processed")
    success_rate: float = Field(..., ge=0.0, le=100.0, description="Percentage of successful operations")


# =============================== Setup Logging ===============================
def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove() # Remove default logger
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )

# =============================== Core Functions ===============================
# DO NOT COPY THESE FUNCTIONS. THEY ARE PLACEHOLDERS.
# CREATE YOUR OWN FUNCTIONS WITH MEANINGFUL NAMES FOR THE TASK.

@logger.catch
def example_placeholder_function(input_data: str) -> float:
    """
    This is an example function. You MUST replace it with a function
    that solves the user's actual problem.
    """
    # This logic is an example. You MUST write your own logic.
    return 42.0

@logger.catch 
def another_placeholder_function(config: InputOutputConfig) -> AnalysisResult:
    """
    This is an example function. You MUST replace it with the core logic
    that solves the user's actual problem.
    """
    # This logic is an example. You MUST write your own logic.
    result = AnalysisResult(
        total_items_processed=100,
        success_rate=95.5
    )
    return result

# =============================== Main Function ===============================
def parse_arguments():
    """Set and parse command line arguments. Modify flags and help text as needed."""
    parser = argparse.ArgumentParser(description="[Brief description of the script's purpose]")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Path to the input file")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

@logger.catch
def main():
    """Main entry point of the script. Replace this docstring with a relevant one."""

    args = parse_arguments()
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)

    # Validate input and output paths using the Pydantic model
    try:
        config = InputOutputConfig(
            input_file=args.input,
            output_file=args.output
        )

        logger.info(f"Starting processing of {config.input_file}")

        # Run the core analysis/logic. REPLACE THIS WITH YOUR LOGIC.
        results = another_placeholder_function(config)
        
        # Save results. MODIFY THIS TO SAVE YOUR SPECIFIC RESULTS.
        if results:
            with open(config.output_file, 'w') as f:
                f.write(results.model_dump_json(indent=2)) # Use .model_dump_json() for Pydantic v2
            logger.success(f"Analysis complete. Results saved to {config.output_file}")
    
    except ValueError as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

```

# Reminder
!!! POST-CODE REVIEW CHECKLIST FOR THE AGENT !!!
1.  **Originality**: Have you replaced ALL placeholder names, logic, and comments with content that directly solves the user’s actual task?
2.  **Docstrings**: 
    *   **If a template/docstring example is given in the prompt**, have you kept **the exact style and format** provided in that template/example?
    *   **If no template/example is given**, have you written **NumPy-style** docstrings for every public API (classes, methods, functions used by external code)?
3.  **Structure**: Are imports, logging setup, Pydantic models, core functions, argument parsing, and the `main()` block in the **same order and format** as the template?
4.  **Error Handling**: Have you consistently used `@logger.catch` (or equivalent) on functions that might raise?
5.  **Logging**: Is every significant step logged at an appropriate level (`info`, `debug`, `success`, `error`)?
6.  **Cleanup**: Have you removed any unused imports, commented-out code, or leftover TODOs?

!!! END CHECKLIST !!!