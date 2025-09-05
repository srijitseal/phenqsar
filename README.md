# PhenQSAR

A phenotypic drug discovery platform for identifying potential therapeutic compounds through Cell Painting analysis.

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://phenqsar.streamlit.app/)
[![Website](https://img.shields.io/badge/Website-srijitseal.com%2Fphenqsar-blue)](https://srijitseal.com/phenqsar/)

## Overview

PhenQSAR finds potential drug candidates by comparing how compounds affect cells. The platform compares compound-induced cellular changes against genetic perturbation profiles (CRISPR knockouts and ORF overexpression) to identify molecules that either mimic or oppose genetic effects.

## Features

- **Gene-Based Search**: Enter a gene name to find compounds with similar or opposite cellular effects
- **Dual Analysis**: Compare against both gene knockout (CRISPR) and overexpression (ORF) data  
- **Chemical Structure Matching**: Identify compounds structurally similar to known drugs
- **Interactive Visualization**: View molecular structures and similarity scores
- **Known Drug Integration**: Cross-reference results with established therapeutic compounds
- **Downloadable Results**: Export analysis data as CSV files

## Demo

**[Try the Live Demo](https://phenqsar.streamlit.app/)**

**[View Documentation](https://srijitseal.github.io/phenqsar/)**

## Quick Start

### Requirements

- Python 3.10+
- Dependencies listed in `requirements.txt`

### Installation

1. Clone the repository:
```bash
git clone https://github.com/srijitseal/phenqsar.git
cd phenqsar/Web_App
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Run the application:
```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`

## Dataset

The demo uses a mini dataset with:
- **40 Cell Painting Features**
- **138K+ Profiles** (compounds, CRISPR, ORF)
- **116K+ Compounds** with metadata
- **11K+ Known Drug Hits**

### Data Setup

Data downloads automatically from Google Drive.

## Usage

1. **Select Gene**: Choose from available genes with CRISPR/ORF data
2. **View Results**:
   - CRISPR knockout hits (mimic/oppose gene deletion effects)
   - ORF overexpression hits (mimic/oppose gene overexpression)
   - Chemical similarity to known drugs
3. **Export**: Download results as CSV

## Analysis Settings

| Parameter | Value |
|-----------|--------|
| Minimum similarity | \|score\| > 0.2 |
| CRISPR results | Top/bottom 100 compounds |
| ORF results | Top/bottom 20 compounds |
| Molecular weight filter | < 1000 Da |
| High similarity threshold | > 0.5 |
| Chemical similarity threshold | > 0.55 |

## Technical Details

### Similarity Scoring
- **Cell Similarity**: Cosine similarity (-1 to +1)
- **Chemical Similarity**: Tanimoto similarity (0 to 1)

### Data Processing Pipeline
The dataset undergoes preprocessing steps indicated in the filename:
- `wellpos` - Position correction
- `var` - Variance adjustment
- `mad` - Outlier detection  
- `outlier` - Outlier removal
- `featselect` - Feature selection
- `sphering` - Data normalization
- `harmony` - Batch correction

### Technology Stack
- **Frontend**: Streamlit
- **Data**: pandas, numpy
- **Chemistry**: RDKit
- **ML**: scikit-learn
- **Visualization**: matplotlib

## Project Structure

```
phenqsar/Web_App/
├── app.py                    # Main application
├── requirements.txt          # Dependencies
├── logo.png                 # App logo
├── docs/                    # GitHub Pages site
└── README.md               # This file
```

## Contributing

1. Fork the repository
2. Create feature branch (`git checkout -b feature/name`)
3. Commit changes (`git commit -am 'Add feature'`)
4. Push branch (`git push origin feature/name`)
5. Open Pull Request

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Support

- **Issues**: Report bugs via [GitHub Issues](https://github.com/srijitseal/phenqsar/issues)
- **Documentation**: Visit the [docs site](https://srijitseal.github.io/phenqsar/)
- **Questions**: Open a discussion on GitHub

## Acknowledgments

- Cell Painting methodology
- Streamlit framework
- RDKit cheminformatics toolkit
- Open source community
