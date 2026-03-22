# PM2.5 Prediction & Uncertainty Quantification (Scotland)

## 🚀 Project Overview
This project develops machine learning models to predict monthly PM2.5 concentrations across Scotland using environmental, meteorological, and spatial data.

The goal is to build a **deployment-ready prediction pipeline** with reliable uncertainty estimation.

---

## 📊 Models Used
- Linear Model (LM)
- LASSO Regression
- Random Forest (RF)
- Quantile Regression Forest (QRF)

---

## 🧪 Validation Strategy
- Leave-One-Site-Out (LOSO)
- Site-level hold-out (80/20 split)

This ensures strong spatial generalisation and avoids overfitting.

---

## 📈 Results
- Best model: **Random Forest / QRF**
- RMSE: **1.43 μg/m³**
- MAE: **0.96 μg/m³**
- R²: **0.658**

---

## 🔍 Key Insights
- Rainfall & wind reduce PM2.5 (removal & dispersion)
- Pressure & humidity increase PM2.5 (stagnation)
- Strong seasonal patterns
- Nonlinear relationships captured by RF

---

## 🧠 Skills Demonstrated
- Data preprocessing & feature engineering
- Machine learning modeling
- Model evaluation (RMSE, MAE, R²)
- Uncertainty quantification (Conformal Prediction)
- Spatial validation

---

## 🛠 Tech Stack
- R / Python
- scikit-learn / ranger
- ggplot2 / matplotlib

---

## 📌 Future Work
- Model deployment pipeline
- Drift detection
- Real-time prediction system
