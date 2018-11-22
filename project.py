import math

# Data

unilateral_familial = {age_one: 4.0, age_two: 6.0, age_three: 8.0,
                       age_four: 4.0, age_five: 2.0, age_older_older: 0}

bilateral_familial = {age_one: 45.0, age_two: 29.0, age_three: 11.0,
                      age_four: 2.0, age_five: 1.0, age_older_older: 2.0}

unilateral_familial = {age_one: 33.0, age_two: 39.0, age_three: 40.0,
                       age_four: 20.0, age_five: 11.0, age_older_older: 6.0}


# Time intervale -> [o,t]
mean_tumors_infinity = 3.0
mean_tumors = -1 * (math.log(math.exp(-1.0 * mean_tumors_infinity) +
                             (1 - math.exp(mean_tumors_infinity)) * herditary_cancer))
probability_eye_no_mutation = float(math.exp(-1 * (float(mean_tumors / 2))))
probability_eye_at_least_one_mutation = 1.0 - probability_no_mutation

probability_individual_no_mutation = float(math.exp(-1.0 * mean_tumors))
probability_individual_at_least_one_mutation = probability_eye_at_least_one_mutation


hereditary_cancer_num = probability_individual_no_mutation - \
    math.exp(-1 * mean_tumors_infinity)
hereditary_cancer_den = 1 - math.exp(-1 * mean_tumors_infinity)
hereditary_cancer = hereditary_cancer_num/float(hereditary_cancer_den)


hereditary_unilat_num = (probability_eye_no_mutation) - \
    (math.exp(-(mean_tumors_infinity/2)))
hereditary_unilat_den = 1 - (math.exp(-(mean_tumors_infinity/2)))

hereditary_unilateral = hereditary_unilat_num/float(hereditary_unilat_den)


hereditary_bilateral_num = (hereditary_unilat_num)**2
hereditary_bilateral_den = (hereditary_unilat_den)**2
hereditary_bilateral = hereditary_bilateral_num/float(hereditary_bilateral_den)
