const std = @import("std");
const Allocator = std.mem.Allocator;

const common = @import("common.zig");
const compute = @import("compute.zig");
const io = @import("io.zig");

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    // Parse command line
    const usr_input = try io.parseArgs(allocator);
    if (usr_input.help) return;

    // Allocate and initialize particules
    var particules = try io.parseInput(allocator, usr_input.input_file);
    defer particules.deinit(allocator);

    // Allocate Potentials
    var u = try common.UpperMatrix.init(allocator, common.N_particules_total);
    defer u.deinit(allocator);

    // Allocate forces
    var forces: common.Forces = common.Forces{};
    defer forces.deinit(allocator);
    try forces.resize(allocator, common.N_particules_total);
    const forces_slice = forces.slice();

    var moments: common.Moments = undefined;
    var vitesses: common.Vitesses = undefined;
    if (usr_input.vitesse_file) |file| {
        vitesses = try io.parseVitesses(allocator, file);

        moments = common.Moments{};
        try moments.resize(allocator, common.N_particules_total);
    } else {
        // Allocate and initialize momentums
        moments = try compute.computeInitMoments(allocator, 0);

        // Allocate speeds
        vitesses = common.Vitesses{};
        try vitesses.resize(allocator, common.N_particules_total);

        // Initialize speeds
        common.vitessesFromMoments(&vitesses, moments);
    }
    defer vitesses.deinit(allocator);
    defer moments.deinit(allocator);

    // Allocate accelerations
    var accelerations = common.Accelerations{};
    defer accelerations.deinit(allocator);
    try accelerations.resize(allocator, common.N_particules_total);

    setZero(&accelerations);

    var fx: f64 = 0.0;
    var fy: f64 = 0.0;
    var fz: f64 = 0.0;

    // Write the first frame with initial state
    try io.writePDB(usr_input.output_file orelse @constCast("test.pdb"), particules, 0, true);

    const delta_t: f64 = @as(f64, @floatFromInt(usr_input.delta_t)) * 1e-3;

    // Compute the potential and forces of the initial system
    compute.potentialV2(particules, &u, &forces);
    var stdout = std.io.getStdOut().writer();

    const csv = try std.fs.cwd().createFile("examen.csv", std.fs.File.CreateFlags{ .read = true, .truncate = true });
    const csv_writer = csv.writer();
    try csv_writer.print("iterations;temperature;Epotentielle;sommeForces;densite;densiteTotale;stddevDensiteTotale;\n", .{});
    var local_sum: f64 = 0;
    var global_sum: f64 = 0;

    for (1..usr_input.total_t) |iter| {

        // Compute the sum of forces
        for (forces_slice.items(.fx), forces_slice.items(.fy), forces_slice.items(.fz)) |x, y, z| {
            fx += x;
            fy += y;
            fz += z;
        }
        const forces_total = fx + fy + fz;
        if (forces_total > -1e-9 and forces_total < 1e-9) {
            try stdout.print("\x1b[33mForces\x1b[0m: {e: >.3} + {e: >.3} + {e: >.3} = \x1b[32m{e: >.3}\x1b[0m\n", .{ fx, fy, fz, fx + fy + fz });
        } else {
            try stdout.print("\x1b[33mForces\x1b[0m: {e: >.3} + {e: >.3} + {e: >.3} = \x1b[31m{e: >.3}\x1b[0m\n", .{ fx, fy, fz, fx + fy + fz });
        }

        // Compute the sum of potentials
        var u_total: f64 = undefined;
        for (u.array) |elem| {
            u_total = u_total + 2 * elem;
        }
        common.momentsFromVitesses(vitesses, &moments);
        const energies = compute.energieCinetique(moments);

        if (energies.temperature < 310 and energies.temperature > 290) {
            try stdout.print("\x1b[33mE Cinétique\x1b[0m:  {d: >8.3}\t\x1b[33mE Potentielle\x1b[0m:  {d: >8.3}  =  \x1b[33mTotal\x1b[0m: {d: >8.3}\t\x1b[33mTempérature\x1b[0m: \x1b[32m{d: >8.1}\x1b[0m K\n", .{ energies.energie, u_total, u_total + energies.energie, energies.temperature });
        } else {
            try stdout.print("\x1b[33mE Cinétique\x1b[0m:  {d: >8.3}\t\x1b[33mE Potentielle\x1b[0m:  {d: >8.3}  =  \x1b[33mTotal\x1b[0m: {d: >8.3}\t\x1b[33mTempérature\x1b[0m: \x1b[31m{d: >8.1}\x1b[0m K\n", .{ energies.energie, u_total, u_total + energies.energie, energies.temperature });
        }

        // Compute one step of velocity verlet
        try compute.velocityVerlet(delta_t, &particules, &vitesses, &accelerations, &forces, &u);

        common.momentsFromVitesses(vitesses, &moments);
        // The center of mass shouldn't move
        compute.correctionMoments(&moments, energies.energie);
        if (iter < 2000 and iter % usr_input.thermostat_step == 0) {
            // Themostat to go closer to the target temperature
            compute.thermostatBeredsen(&moments, energies.temperature);
        }
        common.vitessesFromMoments(&vitesses, moments);

        // Output the pdb
        if (iter % usr_input.save_step == 0) {
            try stdout.print("Ecriture de l'itération {d}\n", .{iter});
            try io.writePDB(usr_input.output_file orelse @constCast("test.pdb"), particules, iter, iter == 0);
        }

        // Output to CSV for gnuplot
        const density = try compute.computeDensity(particules);
        try csv_writer.print("{d};" ** 7, .{ iter, energies.temperature, u_total, forces_total, density.local[0], density.mean, density.mean_stddev });
        try csv_writer.print("\n", .{});
        local_sum += density.local[0];
        global_sum += density.mean;
    }
    csv.close();

    try stdout.print("Moyenne du sous domaine choisi: {d}\n", .{local_sum / @as(f64, @floatFromInt(usr_input.total_t))});
    try stdout.print("Moyenne des moyennes: {d}\n", .{global_sum / @as(f64, @floatFromInt(usr_input.total_t))});
}

// Set to zero all the elements of an positions array
fn setZero(array: *common.Positions) void {
    const array_slice = array.slice();
    const array_x = array_slice.items(.x);
    const array_y = array_slice.items(.y);
    const array_z = array_slice.items(.z);

    for (array_x, array_y, array_z) |*x, *y, *z| {
        x.* = 0;
        y.* = 0;
        z.* = 0;
    }
}

test "Conversions vitesses moments" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var vitesses = common.Vitesses{};
    defer vitesses.deinit(allocator);
    try vitesses.resize(allocator, common.N_particules_total);

    var moments = try compute.computeInitMoments(allocator, 1);
    defer moments.deinit(allocator);

    var moments2 = common.Moments{};
    defer moments2.deinit(allocator);
    try moments2.resize(allocator, common.N_particules_total);

    common.vitessesFromMoments(&vitesses, moments);
    common.momentsFromVitesses(vitesses, &moments2);

    for (moments.items(.fx), moments2.items(.fx)) |x, x2| {
        try std.testing.expectApproxEqRel(x, x2, 1e-15);
    }
    for (moments.items(.fy), moments2.items(.fy)) |y, y2| {
        try std.testing.expectApproxEqRel(y, y2, 1e-15);
    }
    for (moments.items(.fz), moments2.items(.fz)) |z, z2| {
        try std.testing.expectApproxEqRel(z, z2, 1e-15);
    }
}

test "Initialisation moments" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var moments = try compute.computeInitMoments(allocator, 1);
    defer moments.deinit(allocator);

    const tmp = try compute.energieCinetique(moments);

    try std.testing.expectApproxEqAbs(300.0, tmp.temperature, 1e-12);
}

test "Initialisation vitesses" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var vitesses = common.Vitesses{};
    defer vitesses.deinit(allocator);
    try vitesses.resize(allocator, common.N_particules_total);

    var moments = try compute.computeInitMoments(allocator, 1);
    defer moments.deinit(allocator);

    common.vitessesFromMoments(&vitesses, moments);

    var tmp: f64 = 0.0;
    for (vitesses.items(.x)) |x| {
        tmp += x;
    }
    try std.testing.expectApproxEqAbs(0, tmp, 1e-12);
    tmp = 0.0;
    for (vitesses.items(.y)) |y| {
        tmp += y;
    }
    try std.testing.expectApproxEqAbs(0, tmp, 1e-12);
    tmp = 0.0;
    for (vitesses.items(.z)) |z| {
        tmp += z;
    }
    try std.testing.expectApproxEqAbs(0, tmp, 1e-12);
}
