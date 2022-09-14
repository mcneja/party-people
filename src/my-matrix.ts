export { vec2, vec3, mat4 };

type vec2 = [number, number];
type vec3 = [number, number, number];
type mat4 =
    [number, number, number, number,
     number, number, number, number,
     number, number, number, number,
     number, number, number, number];

namespace vec2 {
    export function create(): vec2 {
        return [0, 0];
    }

    export function clone(v: vec2): vec2 {
        return [v[0], v[1]];
    }

    export function fromValues(x0: number, x1: number): vec2 {
        return [x0, x1];
    }

    export function copy(result: vec2, v: vec2) {
        result[0] = v[0];
        result[1] = v[1];
    }

    export function set(result: vec2, x0: number, x1: number) {
        result[0] = x0;
        result[1] = x1;
    }

    export function add(result: vec2, a: vec2, b: vec2) {
        result[0] = a[0] + b[0];
        result[1] = a[1] + b[1];
    }

    export function subtract(result: vec2, a: vec2, b: vec2) {
        result[0] = a[0] - b[0];
        result[1] = a[1] - b[1];
    }

    export function multiply(result: vec2, a: vec2, b: vec2) {
        result[0] = a[0] * b[0];
        result[1] = a[1] * b[1];
    }

    export function scale(result: vec2, a: vec2, scale: number) {
        result[0] = a[0] * scale;
        result[1] = a[1] * scale;
    }

    export function scaleAndAdd(result: vec2, a: vec2, b: vec2, scale: number) {
        result[0] = a[0] + b[0] * scale;
        result[1] = a[1] + b[1] * scale;
    }

    export function distance(a: vec2, b: vec2): number {
        const x = a[0] - b[0];
        const y = a[1] - b[1];
        return Math.hypot(x, y);
    }

    export function squaredDistance(a: vec2, b: vec2): number {
        const x = a[0] - b[0];
        const y = a[1] - b[1];
        return x * x + y * y;
    }

    export function length(a: vec2): number {
        return Math.hypot(a[0], a[1]);
    }

    export function squaredLength(a: vec2): number {
        const x = a[0];
        const y = a[1];
        return x * x + y * y;
    }

    export function negate(result: vec2, a: vec2) {
        result[0] = -a[0];
        result[1] = -a[1];
    }

    export function dot(a: vec2, b: vec2): number {
        return a[0] * b[0] + a[1] * b[1];
    }

    export function lerp(result: vec2, a: vec2, b: vec2, t: number) {
        result[0] = a[0] + t * (b[0] - a[0]);
        result[1] = a[1] + t * (b[1] - a[1]);
    }

    export function zero(result: vec2) {
        result[0] = 0;
        result[1] = 0;
    }
}

namespace vec3 {
    export function create(): vec3 {
        return [0, 0, 0];
    }

    export function clone(v: vec3): vec3 {
        return [v[0], v[1], v[2]];
    }

    export function fromValues(x0: number, x1: number, x2: number): vec3 {
        return [x0, x1, x2];
    }

    export function copy(result: vec3, v: vec3) {
        result[0] = v[0];
        result[1] = v[1];
        result[2] = v[2];
    }

    export function set(result: vec3, x0: number, x1: number, x2: number) {
        result[0] = x0;
        result[1] = x1;
        result[2] = x2;
    }

    export function add(result: vec3, a: vec3, b: vec3) {
        result[0] = a[0] + b[0];
        result[1] = a[1] + b[1];
        result[2] = a[2] + b[2];
    }

    export function subtract(result: vec3, a: vec3, b: vec3) {
        result[0] = a[0] - b[0];
        result[1] = a[1] - b[1];
        result[2] = a[2] - b[2];
    }

    export function multiply(result: vec3, a: vec3, b: vec3) {
        result[0] = a[0] * b[0];
        result[1] = a[1] * b[1];
        result[2] = a[2] * b[2];
    }

    export function scale(result: vec3, a: vec3, scale: number) {
        result[0] = a[0] * scale;
        result[1] = a[1] * scale;
        result[2] = a[2] * scale;
    }

    export function scaleAndAdd(result: vec3, a: vec3, b: vec3, scale: number) {
        result[0] = a[0] + b[0] * scale;
        result[1] = a[1] + b[1] * scale;
        result[2] = a[2] + b[2] * scale;
    }

    export function distance(a: vec3, b: vec3): number {
        const x = a[0] - b[0];
        const y = a[1] - b[1];
        const z = a[2] - b[2];
        return Math.sqrt(x*x + y*y + z*z);
    }

    export function squaredDistance(a: vec3, b: vec3): number {
        const x = a[0] - b[0];
        const y = a[1] - b[1];
        const z = a[2] - b[2];
        return x * x + y * y + z * z;
    }

    export function length(a: vec3): number {
        return Math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    }

    export function squaredLength(a: vec3): number {
        const x = a[0];
        const y = a[1];
        const z = a[2];
        return x * x + y * y + z * z;
    }

    export function negate(result: vec3, a: vec3) {
        result[0] = -a[0];
        result[1] = -a[1];
        result[2] = -a[2];
    }

    export function dot(a: vec3, b: vec3): number {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    export function lerp(result: vec3, a: vec3, b: vec3, t: number) {
        result[0] = a[0] + t * (b[0] - a[0]);
        result[1] = a[1] + t * (b[1] - a[1]);
        result[2] = a[2] + t * (b[2] - a[2]);
    }

    export function zero(result: vec3) {
        result[0] = 0;
        result[1] = 0;
        result[2] = 0;
    }
}

namespace mat4 {
    export function create(): mat4 {
        return [
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
        ];
    }

    export function copy(result: mat4, a: mat4) {
        result[0] = a[0];
        result[1] = a[1];
        result[2] = a[2];
        result[3] = a[3];
        result[4] = a[4];
        result[5] = a[5];
        result[6] = a[6];
        result[7] = a[7];
        result[8] = a[8];
        result[9] = a[9];
        result[10] = a[10];
        result[11] = a[11];
        result[12] = a[12];
        result[13] = a[13];
        result[14] = a[14];
        result[15] = a[15];
    }

    export function identity(m: mat4) {
        m.fill(0);
        m[0] = 1;
        m[5] = 1;
        m[10] = 1;
        m[15] = 1;
    }

    export function translate(m: mat4, v: vec3) {
        // Pre-multiply m by translation matrix v
        m[0] += m[3] * v[0];
        m[1] += m[3] * v[1];
        m[2] += m[3] * v[2];

        m[4] += m[7] * v[0];
        m[5] += m[7] * v[1];
        m[6] += m[7] * v[2];

        m[8] += m[11] * v[0];
        m[9] += m[11] * v[1];
        m[10] += m[11] * v[2];

        m[12] += m[15] * v[0];
        m[13] += m[15] * v[1];
        m[14] += m[15] * v[2];
    }

    export function scale(m: mat4, s: vec3) {
        // Pre-multiply by scale matrix
        m[0] *= s[0];
        m[1] *= s[1];
        m[2] *= s[2];

        m[4] *= s[0];
        m[5] *= s[1];
        m[6] *= s[2];

        m[8] *= s[0];
        m[9] *= s[1];
        m[10] *= s[2];

        m[12] *= s[0];
        m[13] *= s[1];
        m[14] *= s[2];
    }

    export function rotateX(m: mat4, angle: number) {
        const s = Math.sin(angle);
        const c = Math.cos(angle);

        const m1 = m[1];
        const m2 = m[2];
        const m5 = m[5];
        const m6 = m[6];
        const m9 = m[9];
        const m10 = m[10];
        const m13 = m[13];
        const m14 = m[14];

        m[1] = m1 * c + m2 * s;
        m[2] = m2 * c - m1 * s;
        m[5] = m5 * c + m6 * s;
        m[6] = m6 * c - m5 * s;
        m[9] = m9 * c + m10 * s;
        m[10] = m10 * c - m9 * s;
        m[13] = m13 * c + m14 * s;
        m[14] = m14 * c - m13 * s;
    }

    export function frustum(m: mat4, left: number, right: number, bottom: number, top: number, near: number, far: number) {
        const x = (2 * near) / (right - left);
        const y = (2 * near) / (top - bottom);
        const a = (right + left) / (right - left);
        const b = (top + bottom) / (top - bottom);
        const c = (near + far) / (near - far);
        const d = (2 * near * far) / (near - far);

        const m0 = m[0];
        const m1 = m[1];
        const m2 = m[2];
        const m3 = m[3];
        const m4 = m[4];
        const m5 = m[5];
        const m6 = m[6];
        const m7 = m[7];
        const m8 = m[8];
        const m9 = m[9];
        const m10 = m[10];
        const m11 = m[11];
        const m12 = m[12];
        const m13 = m[13];
        const m14 = m[14];
        const m15 = m[15];

        m[0] = m0 * x + m2 * a;
        m[1] = m1 * y + m2 * b;
        m[2] = m2 * c + m3 * d;
        m[3] = -m2;
        m[4] = m4 * x + m6 * a;
        m[5] = m5 * y + m6 * b;
        m[6] = m6 * c + m7 * d;
        m[7] = -m6;
        m[8] = m8 * x + m10 * a;
        m[9] = m9 * y + m10 * b;
        m[10] = m10 * c + m11 * d;
        m[11] = -m10;
        m[12] = m12 * x + m14 * a;
        m[13] = m13 * y + m14 * b;
        m[14] = m14 * c + m15 * d;
        m[15] = -m14;
    }

    export function ortho(result: mat4, left: number, right: number, bottom: number, top: number, near: number, far: number) {
        const lr = 1 / (left - right);
        const bt = 1 / (bottom - top);
        const nf = 1 / (near - far);
        result[0] = -2 * lr;
        result[1] = 0;
        result[2] = 0;
        result[3] = 0;
        result[4] = 0;
        result[5] = -2 * bt;
        result[6] = 0;
        result[7] = 0;
        result[8] = 0;
        result[9] = 0;
        result[10] = 2 * nf;
        result[11] = 0;
        result[12] = (left + right) * lr;
        result[13] = (top + bottom) * bt;
        result[14] = (far + near) * nf;
        result[15] = 1;
    }
}
