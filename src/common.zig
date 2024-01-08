const std = @import("std");

pub const SimuError = error{
    WrongNumberOfArgs,
    WrongNumberOfFields,
};

pub const N_particules_total = 1000;
pub const N_particules_pair = (N_particules_total * (N_particules_total - 1)) / 2;

pub const Particule = struct {
    x: f64,
    y: f64,
    z: f64,
    mass: i8,
};

pub const Force = struct {
    fx: f64 = 0.0,
    fy: f64 = 0.0,
    fz: f64 = 0.0,
};

pub const Particules = std.MultiArrayList(Particule);
pub const Forces = std.MultiArrayList(Force);

pub const CSR = struct {
    array: []f64,
    col_idx: []usize,
    row_ptr: []usize,

    pub fn init(gpa: std.mem.Allocator, n: usize) !CSR {
        const nb_elem: usize = (n * (n - 1)) / 2;
        const self: CSR = CSR{
            .array = try gpa.alloc(f64, nb_elem),
            .col_idx = try gpa.alloc(usize, nb_elem),
            .row_ptr = try gpa.alloc(usize, n),
        };

        @memset(self.array, 0);

        for (0..n) |i| {
            for (i + 1..n) |j| {
                self.col_idx[nb_elem - (((n - i) * (n - i - 1)) / 2) + (j - i - 1)] = j;
            }
        }

        for (self.row_ptr, 0..) |*row, row_idx| {
            row.* = nb_elem - ((n - row_idx) * ((n - row_idx) - 1)) / 2;
        }

        return self;
    }

    pub fn deinit(self: *CSR, gpa: std.mem.Allocator) void {
        gpa.free(self.array);
        gpa.free(self.col_idx);
        gpa.free(self.row_ptr);
    }

    pub fn get(self: CSR, row: usize, col: usize) f64 {
        if (row < col) {
            const i = self.row_ptr[row];
            std.debug.assert(self.col_idx[i + (col - row - 1)] == col);
            return self.array[i + (col - row - 1)];
        } else if (row > col) {
            const i = self.row_ptr[col];
            std.debug.assert(self.col_idx[i + (row - col - 1)] == row);
            return self.array[i + (row - col - 1)];
        } else unreachable;
    }

    pub fn set(self: *CSR, row: usize, col: usize, value: f64) void {
        if (row < col) {
            const i = self.row_ptr[row];
            std.debug.assert(self.col_idx[i + (col - row - 1)] == col);
            self.array[i + (col - row - 1)] = value;
        } else if (row > col) {
            const i = self.row_ptr[col];
            std.debug.assert(self.col_idx[i + (row - col - 1)] == row);
            self.array[i + (row - col - 1)] = value;
        } else unreachable;
    }

    pub fn print(self: CSR) !void {
        const writer = std.io.getStdOut().writer();

        for (0..N_particules_total) |row| {
            for (0..N_particules_total) |col| {
                if (row != col) {
                    try writer.print("{d: >6.5} ", .{self.get(row, col)});
                } else {
                    try writer.print("{d: >6.5} ", .{0.0});
                }
            }
            _ = try writer.write("\n");
        }
    }
};

pub const N_sym: u8 = 1;
pub const R_cut: f64 = 10.0;

pub const boxDim: f64 = 30.00;
