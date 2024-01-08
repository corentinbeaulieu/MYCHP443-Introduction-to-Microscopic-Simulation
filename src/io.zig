const std = @import("std");
const Allocator = std.mem.Allocator;

const common = @import("common.zig");

pub fn parseInput(allocator: Allocator, file: []u8) !common.Particules {
    const open_file = try std.fs.cwd().openFile(file, std.fs.File.OpenFlags{});
    defer open_file.close();

    const file_size = (try open_file.stat()).size;

    const buf: []u8 = try allocator.alloc(u8, file_size);
    defer allocator.free(buf);

    const reader = open_file.reader();
    _ = try reader.readAll(buf);

    var line_split = std.mem.splitScalar(u8, buf, '\n');
    _ = line_split.first();

    var particules: common.Particules = common.Particules{};
    try particules.ensureTotalCapacity(allocator, common.N_particules_total);

    while (line_split.next()) |line| {
        if (line.len == 0) break;
        var field_split = std.mem.splitScalar(u8, line, ' ');
        var field_number: u4 = 0;
        var particule: common.Particule = undefined;
        while (field_split.next()) |field| {
            if (field.len > 0) {
                switch (field_number) {
                    0 => particule.mass = try std.fmt.parseInt(i8, field, 10),
                    1 => particule.x = try std.fmt.parseFloat(f64, field),
                    2 => particule.y = try std.fmt.parseFloat(f64, field),
                    3 => particule.z = try std.fmt.parseFloat(f64, field),
                    else => return common.SimuError.WrongNumberOfFields,
                }
                field_number = field_number + 1;
            }
        }
        particules.appendAssumeCapacity(particule);
    }

    return particules;
}

pub fn writePDB(file: []u8, particules: common.Particules, iteration: usize, truncate: bool) !void {
    const open_file = try std.fs.cwd().createFile(file, std.fs.File.CreateFlags{ .read = true, .truncate = truncate });
    defer open_file.close();
    try open_file.seekFromEnd(0);

    var buffer: [128]u8 = undefined;
    const header = try std.fmt.bufPrint(&buffer, "CRYST1 {0d: >8} {0d: >8} {0d: >8}  90.00  90.00  90.00 P             1\nMODEL  {1d: >7}\n", .{ common.boxDim, iteration });
    _ = try open_file.write(header);

    const particules_slice = particules.slice();
    for (0..common.N_particules_total) |i| {
        const particule = particules_slice.get(i);
        const ligne = try std.fmt.bufPrint(&buffer, "ATOM  {d: >5}  C           0    {d: >8.3}{d: >8.3}{d: >8.3}                  MRES\n", .{ i + 1, particule.x, particule.y, particule.z });
        _ = try open_file.write(ligne);
    }

    const footer = "TER\nENDMDL\n";
    _ = try open_file.write(footer);
}
