const std = @import("std");

pub const N_DOMAINE: usize = 6;
pub const TOTAL_DOMAINES: usize = N_DOMAINE * N_DOMAINE * N_DOMAINE;

pub const Densite = struct {
    local: [TOTAL_DOMAINES]f64 = [_]f64{0} ** TOTAL_DOMAINES,
    mean: f64 = 0,
    mean_stddev: f64 = 0,

    pub fn set_local(self: *Densite, x: usize, y: usize, z: usize, value: f64) void {
        self.local[x * N_DOMAINE * N_DOMAINE + y * N_DOMAINE + z] = value;
    }
};

pub const SimuError = error{
    WrongNumberOfArgs,
    WrongNumberOfFields,
};

pub const N_particules_total: usize = 1000;
pub const N_particules_pair: usize = (N_particules_total * (N_particules_total - 1)) / 2;
pub const N_dl: f64 = 3 * N_particules_total - 3;

pub const r_star: f64 = 3.0;
pub const epsilon_star: f64 = 0.2;

pub const N_sym: u8 = 27;
pub const R_cut: f64 = 10.0;

pub const boxDim: f64 = 32.00;

pub const mass: f64 = 18.0;
pub const CONSTANTE_R: f64 = 0.00199;

pub const CONVERSION_FORCE: f64 = 0.0001 * 4.186;

pub const T_0: f64 = 300.0;
pub const RAPPORT: f64 = (N_dl * CONSTANTE_R * T_0);
pub var rapport: f64 = undefined;
pub const gamma: f64 = 0.01;

pub const epsilon4 = 2 * epsilon_star;
pub const epsilon48 = 48 * epsilon_star;

const dist2_r_cut = (r_star * r_star) / (R_cut * R_cut);
const dist6_r_cut = dist2_r_cut * dist2_r_cut * dist2_r_cut;
const dist8_r_cut = dist2_r_cut * dist2_r_cut * dist2_r_cut * dist2_r_cut;
const dist12_r_cut = dist6_r_cut * dist6_r_cut;
const dist14_r_cut = dist6_r_cut * dist8_r_cut;
pub const potential_r_cut = epsilon4 * (dist12_r_cut - 2 * dist6_r_cut);
pub const force_r_cut = epsilon48 * (dist14_r_cut - dist8_r_cut) * R_cut;

pub const Position = struct {
    x: f64,
    y: f64,
    z: f64,
};
pub const Vitesse = Position;
pub const Acceleration = Position;

pub const Force = struct {
    fx: f64 = 0.0,
    fy: f64 = 0.0,
    fz: f64 = 0.0,
};
pub const Moment = Force;

pub const Positions = std.MultiArrayList(Position);
pub const Forces = std.MultiArrayList(Force);
pub const Vitesses = std.MultiArrayList(Vitesse);
pub const Accelerations = std.MultiArrayList(Acceleration);
pub const Moments = std.MultiArrayList(Moment);

// Structure to represent a symmetric matrix with a null diagonal
pub const UpperMatrix = struct {
    array: []f64,
    size: usize,

    pub fn init(gpa: std.mem.Allocator, n: usize) !UpperMatrix {
        const nb_elem: usize = (n * (n - 1)) / 2;
        const self: UpperMatrix = UpperMatrix{
            .array = try gpa.alloc(f64, nb_elem),
            .size = n,
        };

        @memset(self.array, 0);

        return self;
    }

    pub fn deinit(self: *UpperMatrix, gpa: std.mem.Allocator) void {
        gpa.free(self.array);
        self.size = 0;
    }

    fn localID(self: UpperMatrix, row: usize, col: usize) usize {
        var elem_per_line = self.size - 1;
        var ret = 0;
        if (row > col) {
            for (0..row) |_| {
                ret += elem_per_line;
                elem_per_line -= 1;
            }
            ret += col - (self.size - (elem_per_line + 1));
        } else if (row < col) {
            for (0..col) |_| {
                ret += elem_per_line;
                elem_per_line -= 1;
            }
            ret += row - (self.size - (elem_per_line + 1));
        } else unreachable;

        return ret;
    }

    pub fn get(self: UpperMatrix, row: usize, col: usize) f64 {
        const ret = self.array[self.localID(row, col)];

        return ret;
    }

    pub fn set(self: *UpperMatrix, row: usize, col: usize, value: f64) void {
        self.array[self.localID(row, col)] = value;
    }

    pub fn print(self: UpperMatrix, writer: std.io.AnyWriter) !void {
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

// Converts the momentums into speeds
pub fn vitessesFromMoments(vitesses: *Vitesses, moments: Moments) void {
    const moments_slice = moments.slice();
    const moments_x = moments_slice.items(.fx);
    const moments_y = moments_slice.items(.fy);
    const moments_z = moments_slice.items(.fz);

    const vitesses_slice = vitesses.slice();
    const vitesses_x = vitesses_slice.items(.x);
    const vitesses_y = vitesses_slice.items(.y);
    const vitesses_z = vitesses_slice.items(.z);

    for (vitesses_x, moments_x) |*v_x, m_x| {
        v_x.* = m_x * mass;
    }
    for (vitesses_y, moments_y) |*v_y, m_y| {
        v_y.* = m_y * mass;
    }
    for (vitesses_z, moments_z) |*v_z, m_z| {
        v_z.* = m_z * mass;
    }
}

// Converts the speeds into momentums
pub fn momentsFromVitesses(vitesses: Vitesses, moments: *Moments) void {
    const moments_slice = moments.slice();
    const moments_x = moments_slice.items(.fx);
    const moments_y = moments_slice.items(.fy);
    const moments_z = moments_slice.items(.fz);

    const vitesses_slice = vitesses.slice();
    const vitesses_x = vitesses_slice.items(.x);
    const vitesses_y = vitesses_slice.items(.y);
    const vitesses_z = vitesses_slice.items(.z);

    for (vitesses_x, moments_x) |v_x, *m_x| {
        m_x.* = v_x / mass;
    }
    for (vitesses_y, moments_y) |v_y, *m_y| {
        m_y.* = v_y / mass;
    }
    for (vitesses_z, moments_z) |v_z, *m_z| {
        m_z.* = v_z / mass;
    }
}
