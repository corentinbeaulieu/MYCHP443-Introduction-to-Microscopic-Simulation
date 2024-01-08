const std = @import("std");
const Allocator = std.mem.Allocator;

const common = @import("common.zig");
const compute = @import("compute.zig");
const io = @import("io.zig");

pub fn main() !void {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const args = try std.process.argsAlloc(allocator);

    if (args.len != 2) {
        std.debug.print("USAGE: {s} input_file\n", .{args[0]});
        return common.SimuError.WrongNumberOfArgs;
    }

    var particules = try io.parseInput(allocator, args[1]);
    defer particules.deinit(allocator);

    var results = try compute.potentialU(allocator, particules);
    defer results.u.deinit(allocator);
    defer results.forces.deinit(allocator);

    // try u.print();

    var u_total: f64 = undefined;
    for (results.u.array) |elem| {
        u_total = u_total + elem;
    }

    var fx: f64 = 0.0;
    var fy: f64 = 0.0;
    var fz: f64 = 0.0;
    const forces = results.forces.slice();

    for (forces.items(.fx), forces.items(.fy), forces.items(.fz)) |x, y, z| {
        fx += x;
        fy += y;
        fz += z;
    }

    std.debug.print("Total potential: {d}\n", .{u_total});
    std.debug.print("Energie: x = {e: >.3},\n\t y = {e: >.3},\n\t z = {e: >.3}\n", .{ fx, fy, fz });

    const output_file: []u8 = @constCast("test.pdb");
    try io.writePDB(allocator, output_file, particules, 100);
}
