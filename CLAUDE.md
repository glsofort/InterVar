# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## About InterVar

InterVar is a bioinformatics tool for clinical interpretation of genetic variants using ACMG-AMP 2015 guidelines. It processes genetic variants and classifies them as Pathogenic, Likely Pathogenic, Uncertain Significance, Likely Benign, or Benign based on 28 evidence codes.

## Core Components

- **Intervar.py** (2390 lines) - Main InterVar tool for variant annotation and ACMG classification
- **modify_InterVar.py** (1100 lines) - Modified InterVar implementation with enhanced functionality
- **modify_InterVar_v2.py** (73 lines) - Newer version of the modification script
- **acmg_classification.py** (303 lines) - ACMG classification logic and rules
- **intervardb/** - Database directory containing reference files (mim2gene.txt, domain files, etc.)

## Running InterVar

### Basic Usage
```bash
# Using config file (recommended)
python Intervar.py --config=config.ini

# Direct command line
python Intervar.py -b hg19 -i input_file --input_type=VCF -o output_prefix
```

### Configuration
- Primary config: `config.ini` - Contains all database paths, ANNOVAR locations, and parameters
- Database location: `intervardb/` directory contains all reference files
- ANNOVAR dependency: Requires ANNOVAR >=2016-02-01 installation

### Input Formats
- **AVinput**: ANNOVAR input format (tab-delimited)
- **VCF**: Single sample VCF files
- **VCF_m**: Multi-sample VCF files

## Testing

### Run Tests
```bash
# Run all tests
python test.py

# Tests are located in tests/ directory
python -m unittest tests.test_modify_InterVar
```

## Development Notes

### Key Dependencies
- Python >=2.6.6 (though code uses modern Python 3 features)
- ANNOVAR for variant annotation
- pandas, cyvcf2 for data processing
- Various bioinformatics databases (OMIM, ClinVar, etc.)

### Architecture
1. **Input Processing**: Handles VCF/AVinput formats, calls ANNOVAR if needed
2. **Evidence Code Evaluation**: Implements 28 ACMG evidence codes (PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7)
3. **Classification**: Applies ACMG rules to assign final pathogenicity classification
4. **Output Generation**: Creates detailed reports with evidence codes and classifications

### Modified Versions
- `modify_InterVar.py` includes additional functionality and evidence restrictions
- Uses ClinVar VCF for validation and comparison
- Implements metrics tracking for classification accuracy
- Supports de novo variant analysis

### Evidence Code System
- **PVS**: Very Strong Pathogenic (PVS1)
- **PS**: Strong Pathogenic (PS1-4) 
- **PM**: Moderate Pathogenic (PM1-6)
- **PP**: Supporting Pathogenic (PP1-5)
- **BA**: Stand-alone Benign (BA1)
- **BS**: Strong Benign (BS1-4)
- **BP**: Supporting Benign (BP1-7)

### Database Files Required
- Reference genome annotation databases
- OMIM files (mim2gene.txt, mim_recessive.txt, etc.)
- Domain databases for PM1 evidence
- Population frequency databases
- ClinVar pathogenic/benign variant lists