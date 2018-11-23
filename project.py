import math

mean_tumors_infinity = 3.0
mean_tumors_infinity_true = 4.0

# Data
# Number of sample size
unilateral_familial = {'age_one': 4.0, 'age_two': 6.0, 'age_three': 8.0,
                       'age_four': 4.0, 'age_five': 2.0, 'age_older': 0}

bilateral_familial = {'age_one': 45.0, 'age_two': 29.0, 'age_three': 11.0,
                      'age_four': 2.0, 'age_five': 1.0, 'age_older': 2.0}

unilateral_non_familial = {'age_one': 33.0, 'age_two': 39.0, 'age_three': 40.0,
                           'age_four': 20.0, 'age_five': 11.0, 'age_older': 6.0}


# Calculated results to be compared with predicted
hc_calc_at_infinity = {'age_one': 1.0, 'age_two': 0.62, 'age_three': 0.33,
                       'age_four': 0.13, 'age_five': 0.051, 'age_older': 0.014}

hc_calc_at_infinity_true = {'age_one': 1.0, 'age_two': 0.58, 'age_three': 0.27,
                            'age_four': 0.10, 'age_five': 0.045, 'age_older': 0.017}

hu_calc_at_infinity = {'age_one': 1.0, 'age_two': 0.74, 'age_three': 0.49,
                       'age_four': 0.24, 'age_five': 0.12, 'age_older': 0.037}

hu_calc_at_infinity_true = {'age_one': 1.0, 'age_two': 0.73, 'age_three': 0.46,
                            'age_four': 0.24, 'age_five': 0.13, 'age_older': 0.037}

hb_calc_at_infinity = {'age_one': 1.0, 'age_two': 0.55, 'age_three': 0.24,
                       'age_four': 0.060, 'age_five': 0.014, 'age_older': 0.0013}

hb_calc_at_infinity_true = {'age_one': 1.0, 'age_two': 0.53, 'age_three': 0.22,
                            'age_four': 0.058, 'age_five': 0.017, 'age_older': 0.0034}


m_t = {'age_one': 0.0, 'age_two': 0.0, 'age_three': 0.0,
       'age_four': 0.0, 'age_five': 0.0, 'age_older': 0.0}

m_t_i = {'age_one': 0.0, 'age_two': 0.0, 'age_three': 0.0,
         'age_four': 0.0, 'age_five': 0.0, 'age_older': 0.0}


for i in hc_calc_at_infinity:
    x = math.exp(-1.0 * float(mean_tumors_infinity))
    y = 1 - x
    z = x + y
    prefinal = z * hc_calc_at_infinity[i]
    final = math.log(prefinal)
    # m_t[i] = (math.log(math.exp(-1.0 * mean_tumors_infinity) +
    #    ((1 - math.exp(mean_tumors_infinity)) * hc_calc_at_infinity[i]))) * -1
    m_t[i] = final

for i in hc_calc_at_infinity_true:
    x = math.exp(-1.0 * float(mean_tumors_infinity_true))
    y = 1 - x
    z = x + y
    prefinal = z * hc_calc_at_infinity_true[i]
    final = math.log(prefinal)
    # m_t[i] = (math.log(math.exp(-1.0 * mean_tumors_infinity) +
    #    ((1 - math.exp(mean_tumors_infinity)) * hc_calc_at_infinity[i]))) * -1
    m_t_i[i] = final

print("When m(t) = 3")
for i in m_t:
    print i, m_t[i]

print("\nWhen m(t) = 4")
for i in m_t_i:
    print i, m_t_i[i]

# Functions

# def m_t(hc):
#     calc_mean = -1 * (math.log(math.exp(-1.0 * mean_tumors_infinity) + (
#         (1 - math.exp(mean_tumors_infinity)) * hc)))
#     return calc_mean


# def m_t_i(hc):
#     calc_mean = -1 * (math.log(math.exp(-1.0 * mean_tumors_infinity_true) + (
#         (1 - math.exp(mean_tumors_infinity_true)) * hc)))
#     return calc_mean


def eye_no_mutation(mean_tumors):
    probability_eye_no_mutation = float(
        math.exp(-1 * (float(mean_tumors / 2))))
    return probability_eye_no_mutation


def eye_at_least_one_mutation(probability_no_mutation):
    probability_eye_at_least_one_mutation = 1.0 - \
        float(probability_no_mutation)
    return probability_eye_at_least_one_mutation


def person_no_mutation(mean_tumors):
    probability_individual_no_mutation = float(math.exp(-1.0 * mean_tumors))
    return probability_individual_no_mutation


def person_at_least_one_mutation(probability_eye_at_least_one_mutation):
    probability_individual_at_least_one_mutation = probability_eye_at_least_one_mutation
    return probability_individual_at_least_one_mutation


def HC(m, mi):
    hereditary_cancer_num = (math.exp(-1.0 * m)) - (math.exp(-1.0 * mi))
    hereditary_cancer_den = 1.0 - (math.exp(-1.0 * mi))
    hereditary_cancer = hereditary_cancer_num/float(hereditary_cancer_den)
    return hereditary_cancer


def HU(m, mi):
    hereditary_unilat_num = (math.exp(-1.0 * (float(m)/2))) - \
        (math.exp(-1.0 * (float(mi)/2)))
    hereditary_unilat_den = 1 - (math.exp(-(float(mi)/2)))
    hereditary_unilateral = hereditary_unilat_num/float(hereditary_unilat_den)
    return hereditary_unilateral


def HB(m, mi):
    hereditary_bilateral_num = (
        (math.exp(-1.0 * (float(m)/2))) - (math.exp(-1.0 * (float(mi)/2))))**2
    hereditary_bilateral_den = (1 - (math.exp(-(float(mi)/2))))**2
    hereditary_bilateral = hereditary_bilateral_num / \
        float(hereditary_bilateral_den)
    return hereditary_bilateral


hc_poisson_infinity = {}
hc_poisson_infinity_true = {}

hu_poisson_infinity = {}
hu_poisson_infinity_true = {}

hb_poisson_infinity = {}
hb_poisson_infinity_true = {}
