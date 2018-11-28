import math
import matplotlib.pyplot as plt

mean_tumors_infinity = 3.0
mean_tumors_infinity_true = 4.0

# Data
# Number of sample size
unilateral_familial = {'1': 4.0, '2': 6.0, '3': 8.0,
                       '4': 4.0, '5': 2.0, '6': 0}

bilateral_familial = {'1': 45.0, '2': 29.0, '3': 11.0,
                      '4': 2.0, '5': 1.0, '6': 2.0}

unilateral_non_familial = {'1': 33.0, '2': 39.0, '3': 40.0,
                           '4': 20.0, '5': 11.0, '6': 6.0}


# Calculated results to be compared with predicted
hc_calc_at_infinity = {'1': 1.0, '2': 0.62, '3': 0.33,
                       '4': 0.13, '5': 0.051, '6': 0.014}

hc_calc_at_infinity_true = {'1': 1.0, '2': 0.58, '3': 0.27,
                            '4': 0.10, '5': 0.045, '6': 0.017}

hu_calc_at_infinity = {'1': 1.0, '2': 0.74, '3': 0.49,
                       '4': 0.24, '5': 0.12, '6': 0.037}

hu_calc_at_infinity_true = {'1': 1.0, '2': 0.73, '3': 0.46,
                            '4': 0.24, '5': 0.13, '6': 0.037}

hb_calc_at_infinity = {'1': 1.0, '2': 0.55, '3': 0.24,
                       '4': 0.060, '5': 0.014, '6': 0.0013}

hb_calc_at_infinity_true = {'1': 1.0, '2': 0.53, '3': 0.22,
                            '4': 0.058, '5': 0.017, '6': 0.0034}


m_t = {'1': 0.0, '2': 0.0, '3': 0.0,
       '4': 0.0, '5': 0.0, '6': 0.0}

m_t_i = {'1': 0.0, '2': 0.0, '3': 0.0,
         '4': 0.0, '5': 0.0, '6': 0.0}


for i in hc_calc_at_infinity:
    x = math.exp(-1.0 * float(mean_tumors_infinity))
    y = 1 - x
    z = x + y
    prefinal = z * hc_calc_at_infinity[i]
    final = math.log(prefinal) * -1.0
    # m_t[i] = (math.log(math.exp(-1.0 * mean_tumors_infinity) +
    #    ((1 - math.exp(mean_tumors_infinity)) * hc_calc_at_infinity[i]))) * -1
    m_t[i] = final

for i in hc_calc_at_infinity_true:
    x = math.exp(-1.0 * float(mean_tumors_infinity_true))
    y = 1 - x
    z = x + y
    prefinal = z * hc_calc_at_infinity_true[i]
    final = math.log(prefinal) * -1.0
    # m_t[i] = (math.log(math.exp(-1.0 * mean_tumors_infinity) +
    #    ((1 - math.exp(mean_tumors_infinity)) * hc_calc_at_infinity[i]))) * -1
    m_t_i[i] = final

# print("When m(t) = 3")
# for i in m_t:
#     print i, m_t[i]

# print("\nWhen m(t) = 4")
# for i in m_t_i:
#     print i, m_t_i[i]


# Functions

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


hu_poisson_infinity = {}
hu_poisson_infinity_true = {}

hb_poisson_infinity = {}
hb_poisson_infinity_true = {}


for i in hu_calc_at_infinity:
    hu_poisson_infinity[i] = HU(m_t[i], mean_tumors_infinity)

for i in hu_calc_at_infinity_true:
    hu_poisson_infinity_true[i] = HU(m_t_i[i], mean_tumors_infinity_true)


for i in hb_calc_at_infinity:
    hb_poisson_infinity[i] = HB(m_t[i], mean_tumors_infinity)

for i in hb_calc_at_infinity_true:
    hb_poisson_infinity_true[i] = HB(m_t_i[i], mean_tumors_infinity_true)


print("When m(t) = 3\nHereditary Unilateral calculated")
for i in hu_calc_at_infinity:
    print i, ": ", hu_poisson_infinity[i]

print("\nWhen m(t) = 4\nHereditary Unilateral calculated")
for i in hu_calc_at_infinity_true:
    print i, ": ", hu_poisson_infinity_true[i]
print("-----------------------------------")
print("\nWhen m(t) = 3\nHereditary Bilateral calculated")
for i in hb_calc_at_infinity:
    print i, ": ", hb_poisson_infinity[i]

print("\nWhen m(t) = 4\nHereditary Bilateral calculated")
for i in hb_calc_at_infinity_true:
    print i, ": ", hb_poisson_infinity_true[i]
