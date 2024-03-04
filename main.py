from sympy.physics.continuum_mechanics.beam import Beam
import sympy
import timeit

def main(retained_height, spacing, pile_diameter, pinned, ka, surcharge, soil_weight, phi, shear, moment):
    # Determine loads onto pile




    # Define pressure equation
    x = sympy.Symbol('x', real=True, positive=True)
    eqn1 = (430 * x + 250) * pile_diameter / 1000 * 0.5
    range1 = 4.5 * pile_diameter / 1000
    eqn2 = 1125 * pile_diameter / 1000 * 0.5


    depth = zero_shear_depth(pile_diameter, shear, x, eqn1, range1, eqn2)

    # Determine additional depth required to form moment couple
    if moment > 0:
        sum_moment = sum_moments_above_zero_shear(x, depth, shear, moment, eqn1, range1,
                                                  eqn2)
        additional_depth = moment_couple(eqn1, x, depth, sum_moment, eqn1, range1, eqn2)
        #return additional_depth + depth
        return sum_moment

    return depth


def moment_couple(eqn, x, depth, sum_moment, *args):
    F = sympy.integrate(eqn, (x, depth, x))
    F = sympy.Eq(F * (x - depth), sum_moment)
    additional_depth = sympy.solve([x > depth, F], x)
    if additional_depth.rhs > args[1] and depth > args[1]:
        # must solve for second equation
        F = sympy.integrate(args[2], (x, args[1], x))
        F = sympy.Eq(F * (x - depth), sum_moment)
        additional_depth = sympy.solve([x > depth, F], x)
        return (additional_depth.rhs - depth) * 2
    elif additional_depth.rhs > args[1] and depth < args[1]:
        # must solve for combination of both
        F = sympy.integrate(args[0], (x, depth, args[1])) + sympy.integrate(args[2], (x, args[1], x))
        F = sympy.Eq(F * (x - depth), sum_moment)
        additional_depth = sympy.solve([x > depth, F], x)
        return (additional_depth.rhs - depth) * 2
    return (additional_depth.rhs - depth) * 2


def sum_moments_above_zero_shear(x, depth, shear, moment, *args):
    # integrate Eqn* x / Eqn
    lower_centroid = depth - integrate2(depth, x, args[0] * x, args[1], args[2] * x) / shear
    # lower_centroid = depth - sympy.integrate(eqn*x, (x, 0 , depth))/shear
    moment_below = lower_centroid * shear
    moment_above = shear * depth + moment
    return -moment_below + moment_above


def zero_shear_depth(pile_diameter, shear, x, *args):
    depth = integrate_function(shear, x, *args)
    depth = [i for i in depth if i >= 0][0]
    return depth


def integrate2(depth, x, *args):
    if depth < args[1]:
        result = sympy.integrate(args[0], (x, 0, depth))
        return result
    else:
        result = sympy.integrate(args[0], (x, 0, args[1]))
        result1 = sympy.integrate(args[2], (x, args[1], depth))
        return result + result1


def integrate_function(equals, x, *args):
    # args = [eqn1, range1, eqn2, etc...]
    result = sympy.solve(sympy.Eq(sympy.integrate(args[0], x), equals), x)
    # TODO recursive function to allow for more equation range arguments
    if len(args) > 1:
        if result[0] > args[1]:
            result = sympy.solve(
                sympy.Eq(sympy.integrate(args[2], x), equals - sympy.integrate(args[0], (x, 0, args[1]))), x)
    return result


def arching_factor(phi):
    # phi in degrees
    f = min(2.5, 0.08 * phi)
    return f


def pa(retained_height, ka, spacing, soil_weight):
    # Force due to retained soil
    pa = 1.25 * retained_height ** 2 / 2 * ka * spacing * soil_weight
    return pa


def pw(retained_height, ka, spacing, surcharge):
    # Force due to surcharge
    pw = 1.5 * surcharge * ka * retained_height * spacing
    return pw


def retained_moment(retained_height, Pa, Pw):
    moment = Pa * retained_height / 3 + Pw * retained_height / 2
    return moment


def retained_shear(Pa, Pw):
    shear = Pa + Pw
    return shear


def pinned_rtw(pinned):
    # Run analysis in terms of L once and then substitute to improve runtime.
    L = sympy.symbols('L', positive = True)
    end = str(L)
    E, I, R1, R2, M2, S, ka = sympy.symbols(f'E, I, R1, R2, M2, S, Ka')
    b = Beam(L, E, I)
    b.apply_load(R2, L, -1)
    b.apply_load(-1.25 * 18 * ka * S, 0, 1, L)
    b.apply_load(-1.5 * 5 * ka * S, 0, 0)
    if pinned == True:
        b.apply_load(R1, 0, -1)
        b.bc_deflection = [(0, 0), (L, 0)]
        start = timeit.default_timer()
        b.solve_for_reaction_loads(R1, R2)
        stop = timeit.default_timer()
        print('Time: ', stop - start)
        p = b.reaction_loads
        print(p)
        return p[R2], 0
    else:
        b.apply_load(M2, L, -2)
        b.solve_for_reaction_loads(R2, M2)
        p = b.reaction_loads
        print(p)
        return p[R2], p[M2]


if __name__ == "__main__":
    # Input Parameters
    start = timeit.default_timer()
    ka = 0.42
    surcharge = 5  # in KPa
    soil_weight = 18  # in KN/m3
    phi = 24  # in degrees
    pinned = False
    shear, moment = pinned_rtw(pinned)
    stop = timeit.default_timer()
    print('Time: ', stop - start)
    Ka, L, S = sorted(shear.free_symbols, key=str)
    for i in range(16,26,1):
        print("spacing", i/10)
        for j in range(20,45,5):
            if pinned:
                x = main(j, i, 600, True, ka, surcharge, soil_weight, phi, shear.subs({Ka: 0.42, L: j / 10, S: i / 10}), 0)
            else:

                x = main(j, i, 600, True, ka, surcharge, soil_weight, phi, shear.subs({Ka: 0.42, L: j / 10, S: i / 10}),
                         moment.subs({Ka: 0.42, L: j / 10, S: i / 10}))
            print("retained height", j/10,"Max moment", x)
    print(moment.subs({Ka: 0.42, L: 4, S: 2}))