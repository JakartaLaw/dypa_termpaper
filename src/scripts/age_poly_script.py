
# Libraries
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures



age = np.array([25,30,35,40,45,50,55,60,65]).reshape(-1,1)

col = np.array([40, 50, 57, 62, 66, 67, 65, 61, 54])
hs = np.array([32, 37, 42, 46, 49, 50, 48, 44, 37])
lhs = np.array([26, 31, 34, 38, 39, 39, 37, 34, 31])

educ = [lhs, hs, col]

# Get predicted income
# pred_income = dict()
# for i in educ:
#     # Fit
#     fit = lm.fit(age_poly, i)
#     for t in range(par.start_age, par.retire_age):
#         age_t = np.array(t).reshape(-1,1)
#         age_t_poly = poly.fit_transform(age_t)
#         pred_income[i].append(lm.predict(age_t_poly))

pred_income = dict()
names = ['<HS', 'HS', 'College']
poly = PolynomialFeatures(3)
age_poly =  poly.fit_transform(age)
lm = LinearRegression()

for i in range(len(educ)):
    # Fit
    fit = lm.fit(age_poly, educ[i])
    pred_income[names[i]] = []
    for t in range(25, 65):
        age_t = np.array(t).reshape(-1,1)
        age_t_poly = poly.fit_transform(age_t)
        pred_t = float(lm.predict(age_t_poly))
        pred_income[names[i]].append(pred_t)
