const std = @import("std");
const Allocator = std.mem.Allocator;

const common = @import("./common.zig");

const r_star: f64 = 3.0;
const epsilon_star: f64 = 0.2;

pub fn potentialU(gpa: Allocator, particules: common.Particules) !struct { u: common.CSR, forces: common.Forces } {
    const u = try common.CSR.init(gpa, common.N_particules_total);

    var distances = try computeSquaredDistances(gpa, particules);
    defer distances.deinit(gpa);

    const epsilon4 = 4 * epsilon_star;
    const epsilon48 = -48 * epsilon_star;

    for (u.array, distances.array) |*elem_u, distance| {
        if (distance < common.R_cut) {
            const dist_2: f64 = (r_star * r_star / distance);
            const dist_6 = dist_2 * dist_2 * dist_2;
            elem_u.* = epsilon4 * (dist_6 * dist_6 - (2 * dist_6));
        } else elem_u.* = 0.0;
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
            if (dist < common.R_cut) {
                const dist_2: f64 = (r_star * r_star / dist);
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

fn computeSquaredDistances(gpa: Allocator, particules: common.Particules) !common.CSR {
    var distances = try common.CSR.init(gpa, common.N_particules_total);

    computeSquaredDifferences(particules.items(.x), &distances);
    computeSquaredDifferences(particules.items(.y), &distances);
    computeSquaredDifferences(particules.items(.z), &distances);

    return distances;
}

fn computeSquaredDifferences(positions: []f64, result: *common.CSR) void {
    for (positions, 0..) |p1, i| {
        for (positions[i + 1 ..], i + 1..) |p2, j| {
            const diff = p1 - p2;
            const previous = result.get(i, j);

            result.set(i, j, (diff * diff) + previous);
        }
    }
}
