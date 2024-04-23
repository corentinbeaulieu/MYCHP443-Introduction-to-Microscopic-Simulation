const std = @import("std");
const Allocator = std.mem.Allocator;

const common = @import("./common.zig");

// Examen Question: Density computation
pub fn computeDensity(positions: common.Positions) !common.Densite {
    const subpart: f64 = common.boxDim / common.N_DOMAINE;
    const volume = subpart * subpart * subpart;
    const base: f64 = -common.boxDim / 2;
    const particules_slice = positions.slice();
    const particules_x = particules_slice.items(.x);
    const particules_y = particules_slice.items(.y);
    const particules_z = particules_slice.items(.z);

    var ret = common.Densite{};

    const boxes_dim = [common.N_DOMAINE]f64{ base, base + subpart, base + 2 * subpart, base + 3 * subpart, base + 4 * subpart, base + 5 * subpart };
    for (boxes_dim, 0..) |x_min, i| {
        const x_max = x_min + subpart;
        for (boxes_dim, 0..) |y_min, j| {
            const y_max = y_min + subpart;
            for (boxes_dim, 0..) |z_min, k| {
                const z_max = z_min + subpart;
                var count: usize = 0;
                for (particules_x, particules_y, particules_z) |x, y, z| {
                    if (x > x_min and x < x_max and y > y_min and y < y_max and z > z_min and z < z_max) {
                        count += 1;
                    }
                }
                ret.set_local(i, j, k, @as(f64, @floatFromInt(count)) / volume);
            }
        }
    }
    var sum: f64 = 0;
    for (ret.local) |local| {
        sum += local;
    }
    ret.mean = sum / common.TOTAL_DOMAINES;
    sum = 0;
    for (ret.local) |local| {
        const tmp = local - ret.mean;
        sum += tmp * tmp;
    }
    ret.mean_stddev = std.math.sqrt(sum / common.TOTAL_DOMAINES);

    return ret;
}

// Initialize the Kinetic momentums at the target temperature
pub fn computeInitMoments(allocator: Allocator, seed: u64) !common.Moments {
    var prng = std.rand.DefaultPrng.init(seed);
    const generator = prng.random();

    var moments = common.Moments{};
    try moments.resize(allocator, common.N_particules_total);

    const moments_slice = moments.slice();
    const moments_x = moments_slice.items(.fx);
    const moments_y = moments_slice.items(.fy);
    const moments_z = moments_slice.items(.fz);

    for (moments_x, moments_y, moments_z) |*x, *y, *z| {
        const c_x = generator.float(f64);
        const c_y = generator.float(f64);
        const c_z = generator.float(f64);

        x.* = if (generator.uintAtMost(u1, 1) == 0) c_x else -c_x;
        y.* = if (generator.uintAtMost(u1, 1) == 0) c_y else -c_y;
        z.* = if (generator.uintAtMost(u1, 1) == 0) c_z else -c_z;
    }

    const energie = energieCinetique(moments);

    const rapport = std.math.sqrt(common.RAPPORT / energie.energie);

    for (moments_x, moments_y, moments_z) |*x, *y, *z| {
        x.* *= rapport;
        y.* *= rapport;
        z.* *= rapport;
    }

    return moments;
}

// Computes the Kinetic Energy and associated temperature from the kinetic momentums
pub fn energieCinetique(moments: common.Moments) struct { energie: f64, temperature: f64 } {
    const prefact = 1 / (2 * common.CONVERSION_FORCE * common.mass);

    const moments_slice = moments.slice();
    const moments_x = moments_slice.items(.fx);
    const moments_y = moments_slice.items(.fy);
    const moments_z = moments_slice.items(.fz);

    var energie: f64 = 0;
    for (moments_x, moments_y, moments_z) |x, y, z| {
        energie += (x * x) + (y * y) + (z * z);
    }

    energie *= prefact;

    const temperature = energie / (common.N_dl * common.CONSTANTE_R);

    return .{ .energie = energie, .temperature = temperature };
}

// Computes one step of Velocity Verlet
pub fn velocityVerlet(delta_t: f64, positions: *common.Positions, vitesses: *common.Vitesses, accelerations: *common.Accelerations, forces: *common.Forces, potentiels: *common.UpperMatrix) !void {
    const dt_div_2 = delta_t / 2;

    const positions_slice = positions.slice();
    const positions_x = positions_slice.items(.x);
    const positions_y = positions_slice.items(.y);
    const positions_z = positions_slice.items(.z);

    const vitesses_slice = vitesses.slice();
    const vitesses_x = vitesses_slice.items(.x);
    const vitesses_y = vitesses_slice.items(.y);
    const vitesses_z = vitesses_slice.items(.z);

    const accelerations_slice = accelerations.slice();
    const accelerations_x = accelerations_slice.items(.x);
    const accelerations_y = accelerations_slice.items(.y);
    const accelerations_z = accelerations_slice.items(.z);

    for (positions_x, vitesses_x, accelerations_x) |*p_x, v_x, a_x| {
        p_x.* += v_x * delta_t + a_x * delta_t * dt_div_2;

        if ((p_x.* < -common.boxDim / 2) or (p_x.* > common.boxDim / 2)) {
            p_x.* = try std.math.mod(f64, p_x.*, common.boxDim);
        }

        if (p_x.* > common.boxDim / 2) {
            p_x.* -= common.boxDim;
        }
    }
    for (positions_y, vitesses_y, accelerations_y) |*p_y, v_y, a_y| {
        p_y.* += v_y * delta_t + a_y * delta_t * dt_div_2;
        if ((p_y.* < -common.boxDim / 2) or (p_y.* > common.boxDim / 2)) {
            p_y.* = try std.math.mod(f64, p_y.*, common.boxDim);
        }

        if (p_y.* > common.boxDim / 2) {
            p_y.* -= common.boxDim;
        }
    }
    for (positions_z, vitesses_z, accelerations_z) |*p_z, v_z, a_z| {
        p_z.* += v_z * delta_t + a_z * delta_t * dt_div_2;
        if ((p_z.* < -common.boxDim / 2) or (p_z.* > common.boxDim / 2)) {
            p_z.* = try std.math.mod(f64, p_z.*, common.boxDim);
        }

        if (p_z.* > common.boxDim / 2) {
            p_z.* -= common.boxDim;
        }
    }

    potentialV2(positions.*, potentiels, forces);

    const forces_slice = forces.slice();
    const forces_x = forces_slice.items(.fx);
    const forces_y = forces_slice.items(.fy);
    const forces_z = forces_slice.items(.fz);

    for (vitesses_x, accelerations_x, forces_x) |*v_x, *a_x, f_x| {
        const old_acc = a_x.*;
        a_x.* = f_x / common.mass;
        v_x.* += (a_x.* + old_acc) * dt_div_2;
    }

    for (vitesses_y, accelerations_y, forces_y) |*v_y, *a_y, f_y| {
        const old_acc = a_y.*;
        a_y.* = f_y / common.mass;
        v_y.* += (a_y.* + old_acc) * dt_div_2;
    }

    for (vitesses_z, accelerations_z, forces_z) |*v_z, *a_z, f_z| {
        const old_acc = a_z.*;
        a_z.* = f_z / common.mass;
        v_z.* += (a_z.* + old_acc) * dt_div_2;
    }
}

// Updates the momentums to recenter the mass center
pub fn correctionMoments(moments: *common.Moments, energie: f64) void {
    _ = energie;
    const moments_slice = moments.slice();
    const moments_x = moments_slice.items(.fx);
    const moments_y = moments_slice.items(.fy);
    const moments_z = moments_slice.items(.fz);

    var p_x: f64 = 0;
    var p_y: f64 = 0;
    var p_z: f64 = 0;

    for (moments_x, moments_y, moments_z) |x, y, z| {
        p_x += x;
        p_y += y;
        p_z += z;
    }

    for (moments_x, moments_y, moments_z) |*x, *y, *z| {
        x.* -= (p_x / common.N_particules_total);
        y.* -= (p_y / common.N_particules_total);
        z.* -= (p_z / common.N_particules_total);
    }
    // const rapport = std.math.sqrt(common.RAPPORT / energie);

    // for (moments_x, moments_y, moments_z) |*x, *y, *z| {
    //     x.* *= rapport;
    //     y.* *= rapport;
    //     z.* *= rapport;
    // }
}

// Updates the momentums using the Beredsen thermostat to help keep the target temperature
pub fn thermostatBeredsen(moments: *common.Moments, temperature: f64) void {
    const moments_slice = moments.slice();
    const moments_x = moments_slice.items(.fx);
    const moments_y = moments_slice.items(.fy);
    const moments_z = moments_slice.items(.fz);

    const thermostat = ((common.T_0 / temperature) - 1) * common.gamma;
    for (moments_x, moments_y, moments_z) |*x, *y, *z| {
        x.* += thermostat * x.*;
        y.* += thermostat * y.*;
        z.* += thermostat * z.*;
    }
}

// First version of the potentials and forces computation
// The idea was to compute all the distances ahead of time to keep a good locality
pub fn potentialV1(gpa: Allocator, particules: common.Positions) !struct { u: common.UpperMatrix, forces: common.Forces } {
    const u = try common.UpperMatrix.init(gpa, common.N_particules_total);

    var distances = try computeSquaredDistances(gpa, particules);
    defer distances.deinit(gpa);

    const epsilon4 = 4 * common.epsilon_star;
    const epsilon48 = -48 * common.epsilon_star;

    for (u.array, distances.array) |*elem_u, distance| {
        if (distance < common.R_cut * common.R_cut) {
            const dist_2: f64 = (common.r_star * common.r_star / distance);
            const dist_6 = dist_2 * dist_2 * dist_2;
            elem_u.* = epsilon4 * (dist_6 * dist_6 - (2 * dist_6));
        }
    }

    var forces: common.Forces = common.Forces{};
    try forces.resize(gpa, common.N_particules_total);
    const forces_slice = forces.slice();
    const forces_x = forces_slice.items(.fx);
    const forces_y = forces_slice.items(.fy);
    const forces_z = forces_slice.items(.fz);

    for (forces_x, forces_y, forces_z) |*fx, *fy, *fz| {
        fx.* = 0.0;
        fy.* = 0.0;
        fz.* = 0.0;
    }

    const particule_slice = particules.slice();
    const particule_x = particule_slice.items(.x);
    const particule_y = particule_slice.items(.y);
    const particule_z = particule_slice.items(.z);

    for (forces_x, forces_y, forces_z, particule_x, particule_y, particule_z, 0..) |*fx, *fy, *fz, x, y, z, i| {
        for (forces_x[i + 1 ..], forces_y[i + 1 ..], forces_z[i + 1 ..], particule_x[i + 1 ..], particule_y[i + 1 ..], particule_z[i + 1 ..], i + 1..) |initial_x, initial_y, initial_z, x_j, y_j, z_j, j| {
            const dist = distances.get(i, j);
            if (dist < common.R_cut * common.R_cut) {
                const dist_2: f64 = (common.r_star * common.r_star / dist);
                const dist_4 = dist_2 * dist_2;
                const dist_6 = dist_2 * dist_2 * dist_2;
                const dist_8 = dist_4 * dist_4;
                const dist_14 = dist_6 * dist_8;
                const partial = epsilon48 * (dist_14 - dist_8);

                var force_j = common.Force{ .fx = initial_x, .fy = initial_y, .fz = initial_z };

                const tmp_x = partial * (x - x_j);
                fx.* += tmp_x;
                force_j.fx -= tmp_x;

                const tmp_y = partial * (y - y_j);
                fy.* += tmp_y;
                force_j.fy -= tmp_y;

                const tmp_z = partial * (z - z_j);
                fz.* += tmp_z;
                force_j.fz -= tmp_z;

                forces.set(j, force_j);
            }
        }
    }

    return .{ .u = u, .forces = forces };
}

// Second and final version of the potentials and forces computation
// The idea is to compute only once the possible pairs of particles
pub fn potentialV2(particules: common.Positions, u: *common.UpperMatrix, forces: *common.Forces) void {
    const particule_slice = particules.slice();

    const forces_slice = forces.slice();
    const forces_x = forces_slice.items(.fx);
    const forces_y = forces_slice.items(.fy);
    const forces_z = forces_slice.items(.fz);

    for (forces_x, forces_y, forces_z) |*fx, *fy, *fz| {
        fx.* = 0.0;
        fy.* = 0.0;
        fz.* = 0.0;
    }

    for (u.array) |*elem| {
        elem.* = 0.0;
    }

    const dim = comptime switch (common.N_sym) {
        1 => [_]f64{0},
        27 => [_]f64{ -common.boxDim, 0, common.boxDim },
        else => @compileError("Number of dimensions not implemented"),
    };

    var i: usize = 0;
    var j: usize = 1;
    for (u.array) |*elem_u| {
        if (j % common.N_particules_total == 0) {
            i += 1;
            j = i + 1;
        }

        const particule_i = particule_slice.get(i);
        const particule_j = particule_slice.get(j);
        for (dim) |dim_x| {
            for (dim) |dim_y| {
                for (dim) |dim_z| {
                    const periodic_particule = common.Position{ .x = particule_j.x + dim_x, .y = particule_j.y + dim_y, .z = particule_j.z + dim_z };

                    const distance = computeSquaredDistance(particule_i, periodic_particule);
                    if (distance < common.R_cut * common.R_cut) {
                        const dist_2: f64 = (common.r_star * common.r_star / distance);
                        const dist_6 = dist_2 * dist_2 * dist_2;
                        const partial = common.epsilon48 * dist_6 * (dist_6 - 1.0) / distance;

                        elem_u.* += common.epsilon4 * dist_6 * (dist_6 - 2.0) - common.potential_r_cut;

                        const tmp_x = partial * (particule_i.x - periodic_particule.x) - common.force_r_cut;
                        const tmp_y = partial * (particule_i.y - periodic_particule.y) - common.force_r_cut;
                        const tmp_z = partial * (particule_i.z - periodic_particule.z) - common.force_r_cut;

                        forces_x[i] += tmp_x;
                        forces_y[i] += tmp_y;
                        forces_z[i] += tmp_z;

                        forces_x[j] -= tmp_x;
                        forces_y[j] -= tmp_y;
                        forces_z[j] -= tmp_z;
                    }
                }
            }
        }
        j += 1;
    }
}

// Third version of the potentials and forces computation
// This version is the naivest implementation possible
// It is meant for debugging purpose
pub fn potentialV3(particules: common.Positions, U: []f64, forces: *common.Forces) void {
    const particule_slice = particules.slice();

    const forces_slice = forces.slice();
    const forces_x = forces_slice.items(.fx);
    const forces_y = forces_slice.items(.fy);
    const forces_z = forces_slice.items(.fz);

    for (forces_x, forces_y, forces_z) |*fx, *fy, *fz| {
        fx.* = 0.0;
        fy.* = 0.0;
        fz.* = 0.0;
    }

    for (U) |*elem| {
        elem.* = 0.0;
    }

    const dim = comptime switch (common.N_sym) {
        1 => [_]f64{0},
        27 => [_]f64{ -common.boxDim, 0, common.boxDim },
        else => @compileError("Number of dimensions not implemented"),
    };

    for (dim) |dim_x| {
        for (dim) |dim_y| {
            for (dim) |dim_z| {
                for (0..common.N_particules_total) |i| {
                    const particule_i = particule_slice.get(i);
                    for (0..common.N_particules_total) |j| {
                        if (i == j) continue;
                        const particule_j = particule_slice.get(j);
                        const periodic_particule = common.Position{ .x = particule_j.x + dim_x, .y = particule_j.y + dim_y, .z = particule_j.z + dim_z };

                        const distance = computeSquaredDistance(particule_i, periodic_particule);
                        if (distance < common.R_cut * common.R_cut) {
                            const dist_2: f64 = (common.r_star * common.r_star / distance);
                            const dist_6 = dist_2 * dist_2 * dist_2;
                            const partial = common.epsilon48 * dist_6 * (dist_6 - 1.0) / distance;

                            U[i * common.N_particules_total + j] += common.epsilon4 * (dist_6 * dist_6 - (2 * dist_6));

                            const tmp_x = partial * (particule_i.x - periodic_particule.x) - common.force_r_cut;
                            const tmp_y = partial * (particule_i.y - periodic_particule.y) - common.force_r_cut;
                            const tmp_z = partial * (particule_i.z - periodic_particule.z) - common.force_r_cut;

                            forces_x[i] += tmp_x;
                            forces_y[i] += tmp_y;
                            forces_z[i] += tmp_z;

                            forces_x[j] -= tmp_x;
                            forces_y[j] -= tmp_y;
                            forces_z[j] -= tmp_z;
                        }
                    }
                }
            }
        }
    }
}

// Computes the squared distance between 2 particules
fn computeSquaredDistance(p1: common.Position, p2: common.Position) f64 {
    const x = p1.x - p2.x;
    const y = p1.y - p2.y;
    const z = p1.z - p2.z;

    return x * x + y * y + z * z;
}

// Computes all the squared distances between all the particules
pub fn computeSquaredDistances(gpa: Allocator, particules: common.Positions) !common.UpperMatrix {
    var distances = try common.UpperMatrix.init(gpa, common.N_particules_total);

    computeSquaredDifferences(particules.items(.x), &distances);
    computeSquaredDifferences(particules.items(.y), &distances);
    computeSquaredDifferences(particules.items(.z), &distances);

    return distances;
}

fn computeSquaredDifferences(positions: []f64, result: *common.UpperMatrix) void {
    for (positions, 0..) |p1, i| {
        for (positions[i + 1 ..], i + 1..) |p2, j| {
            const diff = p1 - p2;
            const previous = result.get(i, j);

            result.set(i, j, (diff * diff) + previous);
        }
    }
}
