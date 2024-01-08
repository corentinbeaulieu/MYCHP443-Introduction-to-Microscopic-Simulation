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

pub fn writePDB(allocator: Allocator, file: []u8, particules: common.Particules, iteration: usize) !void {
    const open_file = std.fs.cwd().openFile(file, std.fs.File.OpenFlags{ .mode = .write_only }) catch |err| blk: {
        switch (err) {
            std.fs.File.OpenError.FileNotFound => break :blk try std.fs.cwd().createFile(file, std.fs.File.CreateFlags{ .read = true }),
            else => return err,
        }
    };
    defer open_file.close();

    const header = try std.fmt.allocPrint(allocator, "CRYST1 {0} {0} {0} 90.00 90.00 90.00 P 1\nMODEL {1}", .{ common.boxDim, iteration });
    defer allocator.free(header);
    _ = try open_file.write(header);

    const ligne = try std.fmt.allocPrint(allocator, "ATOM {d: >4} C 0 {d: >16} {d: >16} {d: >16}", .{ 100, particules.get(1).x, particules.get(1).y, particules.get(1).z });
    defer allocator.free(ligne);

    std.debug.print("{}\n", .{ligne.len});

    const length = try open_file.write(ligne);
    std.debug.print("{}\n", .{length});
}
