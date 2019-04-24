import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

def create_age_poly_dict(education_lvl):
    '''Returns a dictionary for the 3 education groups with income profiles (in $1,000) for ages between
    25 and 65, where the first element of the list corresponds to the income at age 25. This data is
    based on the graph on p. 438'''

    # Load data
    age = np.array([25,30,35,40,45,50,55,60,65]).reshape(-1,1)
    col = np.array([40, 50, 57, 62, 66, 67, 65, 61, 54])
    hs = np.array([32, 37, 42, 46, 49, 50, 48, 44, 37])
    lhs = np.array([26, 31, 34, 38, 39, 39, 37, 34, 31])
    educ = [lhs, hs, col]

    # Setup degree 3 regression and output dict format
    lm = LinearRegression()
    poly = PolynomialFeatures(3)
    age_poly =  poly.fit_transform(age)

    names = ['<HS', 'HS', 'College']
    pred_income = dict()

    # Regression and prediction
    for i in range(len(educ)):
        # Fit
        fit = lm.fit(age_poly, educ[i])
        pred_income[names[i]] = []
        for t in range(25, 65):
            age_t = np.array(t).reshape(-1,1)
            age_t_poly = poly.fit_transform(age_t)
            pred_t = float(lm.predict(age_t_poly))
            pred_income[names[i]].append(pred_t)

    return pred_income[education_lvl]
