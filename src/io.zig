const std = @import("std");
const Allocator = std.mem.Allocator;

const clap = @import("clap");

const common = @import("common.zig");

pub const CliArgs = struct {
    input_file: []const u8,
    vitesse_file: ?[]const u8 = null,
    output_file: ?[]const u8 = null,
    delta_t: usize = 1,
    total_t: usize = 10,
    save_step: usize = 10,
    thermostat_step: usize = 1,
    help: bool = false,
};

// Parse the command line options and arguments
pub fn parseArgs(allocator: Allocator) !CliArgs {
    const params = comptime clap.parseParamsComptime(
        \\-h, --help                        Display this help and exit
        \\-d, --delta-t <TIME>              Choose delta t value (default: 1 fs)
        \\-t, --total-iterations <TIME>     Number of iterations to do (default: 100)
        \\-s, --save-step <N>               Save positions every N step (default: 10)
        \\-m, --thermostat-step <N>         Frequency to update the temperature using Beredsen thermostat (default: 1)
        \\<PATH>                            Path to input then output file
    );
    var diag = clap.Diagnostic{};
    const clap_parsers = comptime .{ .N = clap.parsers.int(usize, 10), .TIME = clap.parsers.int(usize, 10), .PATH = clap.parsers.string };
    var res = clap.parse(clap.Help, &params, clap_parsers, .{
        .diagnostic = &diag,
        .allocator = allocator,
    }) catch |err| {
        diag.report(std.io.getStdErr().writer(), err) catch {};
        clap.help(std.io.getStdErr().writer(), clap.Help, &params, .{}) catch {};
        return err;
    };
    defer res.deinit();

    if (res.args.help != 0) {
        clap.help(std.io.getStdOut().writer(), clap.Help, &params, .{}) catch {};
        return CliArgs{ .input_file = "", .help = true };
    }
    if (res.positionals.len == 0) {
        clap.help(std.io.getStdErr().writer(), clap.Help, &params, .{}) catch {};
        return common.SimuError.WrongNumberOfArgs;
    }

    var cli_args = CliArgs{ .input_file = res.positionals[0] };
    if (res.args.@"delta-t") |dt| {
        cli_args.delta_t = dt;
    }
    if (res.args.@"total-iterations") |total_t| {
        cli_args.total_t = total_t;
    }
    if (res.args.@"save-step") |save_step| {
        cli_args.save_step = save_step;
    }
    if (res.args.@"thermostat-step") |thermostat_step| {
        cli_args.thermostat_step = thermostat_step;
    }
    if (res.positionals.len >= 2) {
        cli_args.vitesse_file = res.positionals[1];
    }
    if (res.positionals.len >= 3) {
        cli_args.output_file = res.positionals[2];
    }

    return cli_args;
}

// Parses the input file with the particules positions
pub fn parseInput(allocator: Allocator, file: []const u8) !common.Positions {
    const open_file = try std.fs.cwd().openFile(file, std.fs.File.OpenFlags{});
    defer open_file.close();

    const file_size = (try open_file.stat()).size;

    const buf: []u8 = try allocator.alloc(u8, file_size);
    defer allocator.free(buf);

    const reader = open_file.reader();
    _ = try reader.readAll(buf);

    var line_split = std.mem.splitScalar(u8, buf, '\n');
    _ = line_split.first();

    var particules: common.Positions = common.Positions{};
    try particules.ensureTotalCapacity(allocator, common.N_particules_total);

    while (line_split.next()) |line| {
        if (line.len == 0) break;
        var field_split = std.mem.splitScalar(u8, line, ' ');
        var field_number: u4 = 0;
        var particule: common.Position = undefined;
        while (field_split.next()) |field| {
            if (field.len > 0) {
                switch (field_number) {
                    0 => _ = field.len, // particule.mass = try std.fmt.parseInt(i8, field, 10),
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

// Parses the input file with the particules speeds
pub fn parseVitesses(allocator: Allocator, file: []const u8) !common.Vitesses {
    const open_file = try std.fs.cwd().openFile(file, std.fs.File.OpenFlags{});
    defer open_file.close();

    const file_size = (try open_file.stat()).size;

    const buf: []u8 = try allocator.alloc(u8, file_size);
    defer allocator.free(buf);

    const reader = open_file.reader();
    _ = try reader.readAll(buf);

    var line_split = std.mem.splitScalar(u8, buf, '\n');
    _ = line_split.first();

    var vitesses: common.Positions = common.Vitesses{};
    try vitesses.ensureTotalCapacity(allocator, common.N_particules_total);

    while (line_split.next()) |line| {
        //if (line.len == 0) break;
        var field_split = std.mem.splitScalar(u8, line, ' ');
        var field_number: u4 = 0;
        var vitesse: common.Vitesse = undefined;
        while (field_split.next()) |field| {
            if (field.len > 0) {
                switch (field_number) {
                    0 => _ = field.len, // particule.mass = try std.fmt.parseInt(i8, field, 10),
                    1 => vitesse.x = try std.fmt.parseFloat(f64, field),
                    2 => vitesse.y = try std.fmt.parseFloat(f64, field),
                    3 => vitesse.z = try std.fmt.parseFloat(f64, field),
                    else => return common.SimuError.WrongNumberOfFields,
                }
                field_number = field_number + 1;
            }
        }
        vitesses.appendAssumeCapacity(vitesse);
    }

    return vitesses;
}

// Store the current positions in the pdb format
pub fn writePDB(file: []const u8, particules: common.Positions, iteration: usize, truncate: bool) !void {
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
