import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path


def walsh(coefficients):
    return [sum((-1 if (frequency & position).bit_count() % 2 else 1) * value
                for position, value in enumerate(coefficients))
            for frequency in range(len(coefficients))]


def log_unit_interval(value, terms=70):
    ratio = (value - 1) / (value + 1)
    partial = sum((2 * ratio ** (2 * index + 1) / (2 * index + 1)
                   for index in range(terms)), Fraction(0))
    tail = 2 * ratio ** (2 * terms + 1) / ((2 * terms + 1) * (1 - ratio * ratio))
    return partial, partial + tail


def log_interval(value):
    value = Fraction(value)
    assert value > 0
    power = 0
    while value < 1:
        value *= 2
        power -= 1
    while value > 2:
        value /= 2
        power += 1
    lower, upper = log_unit_interval(value)
    log_two_lower, log_two_upper = log_unit_interval(Fraction(2))
    if power >= 0:
        return lower + power * log_two_lower, upper + power * log_two_upper
    return lower + power * log_two_upper, upper + power * log_two_lower


def entropy_interval(probabilities):
    lower = Fraction(0)
    upper = Fraction(0)
    for probability in probabilities:
        if probability:
            log_lower, log_upper = log_interval(probability)
            lower -= probability * log_upper
            upper -= probability * log_lower
    return lower, upper


def rational_strings(values):
    return [str(value) for value in values]


def basis_distribution(coefficients, quadratic, linear):
    dimension = len(linear)
    assert len(coefficients) == 1 << dimension
    assert len(quadratic) == dimension * (dimension - 1) // 2
    norm_squared = sum(value * value for value in coefficients)
    numerators = []
    for frequency in range(len(coefficients)):
        real_part = 0
        imaginary_part = 0
        for position, value in enumerate(coefficients):
            coordinates = [(position >> shift) & 1 for shift in reversed(range(dimension))]
            quadratic_exponent = sum(value * coordinates[left] * coordinates[right]
                                     for value, (left, right) in
                                     zip(quadratic, combinations(range(dimension), 2)))
            linear_exponent = sum(left * right for left, right in zip(linear, coordinates))
            sign = -1 if (quadratic_exponent + (frequency & position).bit_count()) % 2 else 1
            phase_real, phase_imaginary = [(1, 0), (0, -1), (-1, 0), (0, 1)][linear_exponent % 4]
            real_part += value * sign * phase_real
            imaginary_part += value * sign * phase_imaginary
        numerators.append(real_part * real_part + imaginary_part * imaginary_part)
    probabilities = [Fraction(value, len(coefficients) * norm_squared) for value in numerators]
    assert sum(probabilities) == 1
    return probabilities


def decimal_text(value):
    with localcontext() as context:
        context.prec = 18
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def check():
    signed = [1, 0, 2, -1, 1, 0, 0, 1]
    absolute = list(map(abs, signed))
    signed_walsh = walsh(signed)
    absolute_walsh = walsh(absolute)
    assert signed_walsh == [4, 4, 0, 0, 0, 4, 0, -4]
    assert absolute_walsh == [6, 2, -2, 2, 2, 2, -2, -2]
    signed_probabilities = basis_distribution(signed, (0, 0, 0), (0, 0, 0))
    absolute_probabilities = basis_distribution(absolute, (0, 0, 0), (0, 0, 0))
    assert signed_probabilities == [Fraction(value * value, 64) for value in signed_walsh]
    assert absolute_probabilities == [Fraction(9, 16)] + [Fraction(1, 16)] * 7
    assert 2 ** 16 > 3 ** 9
    assert sum(value * value for value in absolute_probabilities) == Fraction(11, 32)
    path_probabilities = basis_distribution(absolute, (1, 0, 1), (0, 0, 0))
    assert path_probabilities == signed_probabilities
    for first_sign, second_sign in product((1, -1), repeat=2):
        transformed = walsh([1, 0, 2, first_sign, 1, 0, 0, second_sign])
        probabilities = [Fraction(value * value, 64) for value in transformed]
        expected = absolute_probabilities if first_sign == second_sign else signed_probabilities
        assert sorted(probabilities) == sorted(expected)
    tensor_coefficients = [left * right for left in signed for right in signed]
    tensor_p = [Fraction(value * value, 4096) for value in walsh(tensor_coefficients)]
    tensor_q = [Fraction(value * value, 4096) for value in walsh(list(map(abs, tensor_coefficients)))]
    assert tensor_p == [left * right for left in signed_probabilities for right in signed_probabilities]
    assert tensor_q == [left * right for left in absolute_probabilities for right in absolute_probabilities]
    base_count = 0
    minimizer_count = 0
    log_four_lower, log_four_upper = log_interval(Fraction(4))
    for quadratic in product((0, 1), repeat=3):
        for linear in product((0, 1), repeat=3):
            probabilities = basis_distribution(absolute, quadratic, linear)
            assert sum(value * value for value in probabilities) <= Fraction(11, 32)
            modified_quadratic = tuple(value ^ (linear[left] * linear[right])
                                       for value, (left, right) in
                                       zip(quadratic, combinations(range(3), 2)))
            reduced_probabilities = basis_distribution(absolute, modified_quadratic, (0, 0, 0))
            if sorted(reduced_probabilities) != sorted(probabilities):
                reduced_lower, reduced_upper = entropy_interval(reduced_probabilities)
                original_lower, original_upper = entropy_interval(probabilities)
                assert reduced_upper < original_lower
            positive = [value for value in probabilities if value]
            if sorted(positive) == [Fraction(1, 4)] * 4:
                minimizer_count += 1
            else:
                entropy_lower, entropy_upper = entropy_interval(probabilities)
                assert entropy_lower > log_four_upper
            base_count += 1
    graph_coefficients = [1, 1, 1, -1]
    graph_x = basis_distribution(graph_coefficients, (0,), (0, 0))
    graph_y = basis_distribution(graph_coefficients, (0,), (1, 1))
    graph_real_reduced = basis_distribution(graph_coefficients, (1,), (0, 0))
    assert graph_x == [Fraction(1, 4)] * 4
    assert graph_y == [Fraction(1, 2), Fraction(0), Fraction(0), Fraction(1, 2)]
    assert graph_real_reduced == [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
    strict_coefficients = [100, 1, 200, 100, 100, 1, 1, 100]
    assert sum(value * value for value in strict_coefficients) == 80003
    strict_absolute = basis_distribution(strict_coefficients, (0, 0, 0), (0, 0, 0))
    strict_path = basis_distribution(strict_coefficients, (1, 0, 1), (0, 0, 0))
    strict_q_lower, strict_q_upper = entropy_interval(strict_absolute)
    strict_p_lower, strict_p_upper = entropy_interval(strict_path)
    gap_lower = strict_q_lower - strict_p_upper
    gap_upper = strict_q_upper - strict_p_lower
    assert Fraction(137, 1000) < gap_lower < gap_upper < Fraction(139, 1000)
    assert strict_q_upper - strict_q_lower < Fraction(1, 10 ** 50)
    assert strict_p_upper - strict_p_lower < Fraction(1, 10 ** 50)
    for dimension in range(3, 9):
        padded = signed + [0] * ((1 << dimension) - 8)
        transformed = walsh(padded)
        norm = sum(value * value for value in padded)
        probabilities = [Fraction(value * value, (1 << dimension) * norm) for value in transformed]
        assert sum(probabilities) == 1
        assert sorted(value for value in probabilities if value) == [Fraction(1, 1 << (dimension - 1))] * (1 << (dimension - 1))
    return {
        'status': 'pass',
        'coefficient_order': '000,001,010,011,100,101,110,111',
        'signed_coefficients': signed,
        'normalization_squared': 8,
        'signed_walsh': signed_walsh,
        'absolute_walsh': absolute_walsh,
        'p': rational_strings(signed_probabilities),
        'q': rational_strings(absolute_probabilities),
        'strict_gap_identity': '(16*ln(2)-9*ln(3))/8 > 0',
        'integer_certificate': {'two_power_16': 2 ** 16, 'three_power_9': 3 ** 9},
        'positive_amplitude_phase': 'x0*x1+x1*x2',
        'fixed_R_real_phase_sign_pairs_checked': 4,
        'tensor_product_copies_checked': 2,
        'fixed_R_identity_enumeration': {
            'quadratic_and_linear_choices': base_count,
            'minimum_entropy': 'ln(4)',
            'minimizers': minimizer_count,
            'method': 'Exact Gaussian-integer overlaps and rational logarithm bounds',
            'scope': 'The nonnegative three-qubit state; R=I3; all 8 quadratic and 8 linear choices'
        },
        'strictly_positive_example': {
            'coefficients': strict_coefficients,
            'normalization_squared': 80003,
            'q': rational_strings(strict_absolute),
            'p': rational_strings(strict_path),
            'certified_gap_lower': '137/1000',
            'certified_gap_upper': '139/1000',
            'illustrative_decimal_gap_lower': decimal_text(gap_lower),
            'illustrative_decimal_gap_upper': decimal_text(gap_upper),
            'entropy_interval_width_bound': '1e-50',
            'logarithm_series_terms': 70
        },
        'integer_sum_fixed_quadratic_counterexample': {
            'coefficients': graph_coefficients,
            'normalization_squared': 4,
            'Qprime': '0',
            'c_zero_distribution': rational_strings(graph_x),
            'c_11_distribution': rational_strings(graph_y),
            'corrected_real_basis_distribution': rational_strings(graph_real_reduced),
            'entropy_comparison': 'ln(4) > ln(2) for every Renyi index alpha>=0',
            'corrected_reduction_shannon_choices_checked_at_three_qubits': base_count
        },
        'padding_finite_checks_dimensions': list(range(3, 9)),
        'written_proofs': ['Shannon gap', 'Nonnegative-amplitude extension',
                           'Strictly positive existence by continuity', 'All dimensions N>=3',
                           'Renyi interval defined by the unique root alpha_star>1',
                           'Extensive Shannon defect',
                           'Real-amplitude reduction with explicit quadratic phase conversion',
                           'All-index fixed-R minimum for the witness',
                           'Collision-entropy lower bound for arbitrary affine supports',
                           'Coset-conditioned Shannon lower bound',
                           'All-index validity on coset-flat amplitudes and Figure 1 inputs'],
        'lean_declaration': 'D5.S3.Quantum.Magic.StabilizerAsymmetryShannonBoundRefutation.result',
        'lean_scope': 'Negation of the full-support Shannon bound, witnessed by the signed three-qubit state'
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Exact finite checks for the stabilizer Shannon bound.')
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    result = json.dumps(check(), indent=2) + '\n'
    if arguments.output:
        arguments.output.write_text(result)
    else:
        print(result, end='')
