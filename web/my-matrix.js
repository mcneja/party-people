export { vec2, vec3, mat4, MatrixStack };
var vec2;
(function (vec2) {
    function create() {
        return [0, 0];
    }
    vec2.create = create;
    function clone(v) {
        return [v[0], v[1]];
    }
    vec2.clone = clone;
    function fromValues(x0, x1) {
        return [x0, x1];
    }
    vec2.fromValues = fromValues;
    function copy(result, v) {
        result[0] = v[0];
        result[1] = v[1];
    }
    vec2.copy = copy;
    function set(result, x0, x1) {
        result[0] = x0;
        result[1] = x1;
    }
    vec2.set = set;
    function add(result, a, b) {
        result[0] = a[0] + b[0];
        result[1] = a[1] + b[1];
    }
    vec2.add = add;
    function subtract(result, a, b) {
        result[0] = a[0] - b[0];
        result[1] = a[1] - b[1];
    }
    vec2.subtract = subtract;
    function multiply(result, a, b) {
        result[0] = a[0] * b[0];
        result[1] = a[1] * b[1];
    }
    vec2.multiply = multiply;
    function scale(result, a, scale) {
        result[0] = a[0] * scale;
        result[1] = a[1] * scale;
    }
    vec2.scale = scale;
    function scaleAndAdd(result, a, b, scale) {
        result[0] = a[0] + b[0] * scale;
        result[1] = a[1] + b[1] * scale;
    }
    vec2.scaleAndAdd = scaleAndAdd;
    function distance(a, b) {
        const x = a[0] - b[0];
        const y = a[1] - b[1];
        return Math.hypot(x, y);
    }
    vec2.distance = distance;
    function squaredDistance(a, b) {
        const x = a[0] - b[0];
        const y = a[1] - b[1];
        return x * x + y * y;
    }
    vec2.squaredDistance = squaredDistance;
    function length(a) {
        return Math.hypot(a[0], a[1]);
    }
    vec2.length = length;
    function squaredLength(a) {
        const x = a[0];
        const y = a[1];
        return x * x + y * y;
    }
    vec2.squaredLength = squaredLength;
    function negate(result, a) {
        result[0] = -a[0];
        result[1] = -a[1];
    }
    vec2.negate = negate;
    function dot(a, b) {
        return a[0] * b[0] + a[1] * b[1];
    }
    vec2.dot = dot;
    function perpDot(a, b) {
        return a[0] * b[1] - a[1] * b[0];
    }
    vec2.perpDot = perpDot;
    function lerp(result, a, b, t) {
        result[0] = a[0] + t * (b[0] - a[0]);
        result[1] = a[1] + t * (b[1] - a[1]);
    }
    vec2.lerp = lerp;
    function zero(result) {
        result[0] = 0;
        result[1] = 0;
    }
    vec2.zero = zero;
})(vec2 || (vec2 = {}));
var vec3;
(function (vec3) {
    function create() {
        return [0, 0, 0];
    }
    vec3.create = create;
    function clone(v) {
        return [v[0], v[1], v[2]];
    }
    vec3.clone = clone;
    function fromValues(x0, x1, x2) {
        return [x0, x1, x2];
    }
    vec3.fromValues = fromValues;
    function copy(result, v) {
        result[0] = v[0];
        result[1] = v[1];
        result[2] = v[2];
    }
    vec3.copy = copy;
    function set(result, x0, x1, x2) {
        result[0] = x0;
        result[1] = x1;
        result[2] = x2;
    }
    vec3.set = set;
    function add(result, a, b) {
        result[0] = a[0] + b[0];
        result[1] = a[1] + b[1];
        result[2] = a[2] + b[2];
    }
    vec3.add = add;
    function subtract(result, a, b) {
        result[0] = a[0] - b[0];
        result[1] = a[1] - b[1];
        result[2] = a[2] - b[2];
    }
    vec3.subtract = subtract;
    function multiply(result, a, b) {
        result[0] = a[0] * b[0];
        result[1] = a[1] * b[1];
        result[2] = a[2] * b[2];
    }
    vec3.multiply = multiply;
    function scale(result, a, scale) {
        result[0] = a[0] * scale;
        result[1] = a[1] * scale;
        result[2] = a[2] * scale;
    }
    vec3.scale = scale;
    function scaleAndAdd(result, a, b, scale) {
        result[0] = a[0] + b[0] * scale;
        result[1] = a[1] + b[1] * scale;
        result[2] = a[2] + b[2] * scale;
    }
    vec3.scaleAndAdd = scaleAndAdd;
    function distance(a, b) {
        const x = a[0] - b[0];
        const y = a[1] - b[1];
        const z = a[2] - b[2];
        return Math.sqrt(x * x + y * y + z * z);
    }
    vec3.distance = distance;
    function squaredDistance(a, b) {
        const x = a[0] - b[0];
        const y = a[1] - b[1];
        const z = a[2] - b[2];
        return x * x + y * y + z * z;
    }
    vec3.squaredDistance = squaredDistance;
    function length(a) {
        return Math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    }
    vec3.length = length;
    function squaredLength(a) {
        const x = a[0];
        const y = a[1];
        const z = a[2];
        return x * x + y * y + z * z;
    }
    vec3.squaredLength = squaredLength;
    function negate(result, a) {
        result[0] = -a[0];
        result[1] = -a[1];
        result[2] = -a[2];
    }
    vec3.negate = negate;
    function dot(a, b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
    vec3.dot = dot;
    function lerp(result, a, b, t) {
        result[0] = a[0] + t * (b[0] - a[0]);
        result[1] = a[1] + t * (b[1] - a[1]);
        result[2] = a[2] + t * (b[2] - a[2]);
    }
    vec3.lerp = lerp;
    function zero(result) {
        result[0] = 0;
        result[1] = 0;
        result[2] = 0;
    }
    vec3.zero = zero;
    function transformPosition(result, a, m) {
        const x = a[0];
        const y = a[1];
        const z = a[2];
        const w = m[3] * x + m[7] * y + m[11] * z + m[15];
        result[0] = (m[0] * x + m[4] * y + m[8] * z + m[12]) / w;
        result[1] = (m[1] * x + m[5] * y + m[9] * z + m[13]) / w;
        result[2] = (m[2] * x + m[6] * y + m[10] * z + m[14]) / w;
    }
    vec3.transformPosition = transformPosition;
    function transformDirection(result, a, m) {
        const x = a[0];
        const y = a[1];
        const z = a[2];
        const w = m[3] * x + m[7] * y + m[11] * z;
        result[0] = (m[0] * x + m[4] * y + m[8] * z) / w;
        result[1] = (m[1] * x + m[5] * y + m[9] * z) / w;
        result[2] = (m[2] * x + m[6] * y + m[10] * z) / w;
    }
    vec3.transformDirection = transformDirection;
})(vec3 || (vec3 = {}));
var mat4;
(function (mat4) {
    function create() {
        return [
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
        ];
    }
    mat4.create = create;
    function clone(a) {
        return [
            a[0], a[1], a[2], a[3],
            a[4], a[5], a[6], a[7],
            a[8], a[9], a[10], a[11],
            a[12], a[13], a[14], a[15],
        ];
    }
    mat4.clone = clone;
    function copy(result, a) {
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
    mat4.copy = copy;
    function identity(m) {
        m.fill(0);
        m[0] = 1;
        m[5] = 1;
        m[10] = 1;
        m[15] = 1;
    }
    mat4.identity = identity;
    function transpose(result, m) {
        const m00 = m[0];
        const m10 = m[1];
        const m20 = m[2];
        const m30 = m[3];
        const m01 = m[4];
        const m11 = m[5];
        const m21 = m[6];
        const m31 = m[7];
        const m02 = m[8];
        const m12 = m[9];
        const m22 = m[10];
        const m32 = m[11];
        const m03 = m[12];
        const m13 = m[13];
        const m23 = m[14];
        const m33 = m[15];
        result[0] = m00;
        result[1] = m01;
        result[2] = m02;
        result[3] = m03;
        result[4] = m10;
        result[5] = m11;
        result[6] = m12;
        result[7] = m13;
        result[8] = m20;
        result[9] = m21;
        result[10] = m22;
        result[11] = m23;
        result[12] = m30;
        result[13] = m31;
        result[14] = m32;
        result[15] = m33;
    }
    mat4.transpose = transpose;
    function translate(m, v) {
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
    mat4.translate = translate;
    function scale(m, s) {
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
    mat4.scale = scale;
    function rotateX(m, angle) {
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
    mat4.rotateX = rotateX;
    function rotateZ(m, angle) {
        const s = Math.sin(angle);
        const c = Math.cos(angle);
        const m0 = m[0];
        const m1 = m[1];
        const m4 = m[4];
        const m5 = m[5];
        const m8 = m[8];
        const m9 = m[9];
        const m12 = m[12];
        const m13 = m[13];
        m[0] = m0 * c - m1 * s;
        m[1] = m0 * s + m1 * c;
        m[4] = m4 * c - m5 * s;
        m[5] = m4 * s + m5 * c;
        m[8] = m8 * c - m9 * s;
        m[9] = m8 * s + m9 * c;
        m[12] = m12 * c - m13 * s;
        m[13] = m12 * s + m13 * c;
    }
    mat4.rotateZ = rotateZ;
    function frustum(m, left, right, bottom, top, near, far) {
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
    mat4.frustum = frustum;
    function ortho(result, left, right, bottom, top, near, far) {
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
    mat4.ortho = ortho;
    function multiply(result, a, b) {
        const aCopy = [...a];
        const bCopy = [...b];
        for (let col = 0; col < 4; ++col) {
            for (let row = 0; row < 4; ++row) {
                let x = 0;
                for (let i = 0; i < 4; ++i) {
                    x += aCopy[row + 4 * i] * bCopy[4 * col + i];
                }
                result[4 * col + row] = x;
            }
        }
    }
    mat4.multiply = multiply;
    function invert(result, a) {
        const a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3];
        const a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7];
        const a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11];
        const a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15];
        const b00 = a00 * a11 - a01 * a10;
        const b01 = a00 * a12 - a02 * a10;
        const b02 = a00 * a13 - a03 * a10;
        const b03 = a01 * a12 - a02 * a11;
        const b04 = a01 * a13 - a03 * a11;
        const b05 = a02 * a13 - a03 * a12;
        const b06 = a20 * a31 - a21 * a30;
        const b07 = a20 * a32 - a22 * a30;
        const b08 = a20 * a33 - a23 * a30;
        const b09 = a21 * a32 - a22 * a31;
        const b10 = a21 * a33 - a23 * a31;
        const b11 = a22 * a33 - a23 * a32;
        let det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
        if (!det) {
            copy(result, a);
            return;
        }
        det = 1.0 / det;
        result[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
        result[1] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
        result[2] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
        result[3] = (a22 * b04 - a21 * b05 - a23 * b03) * det;
        result[4] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
        result[5] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
        result[6] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
        result[7] = (a20 * b05 - a22 * b02 + a23 * b01) * det;
        result[8] = (a10 * b10 - a11 * b08 + a13 * b06) * det;
        result[9] = (a01 * b08 - a00 * b10 - a03 * b06) * det;
        result[10] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
        result[11] = (a21 * b02 - a20 * b04 - a23 * b00) * det;
        result[12] = (a11 * b07 - a10 * b09 - a12 * b06) * det;
        result[13] = (a00 * b09 - a01 * b07 + a02 * b06) * det;
        result[14] = (a31 * b01 - a30 * b03 - a32 * b00) * det;
        result[15] = (a20 * b03 - a21 * b01 + a22 * b00) * det;
    }
    mat4.invert = invert;
})(mat4 || (mat4 = {}));
class MatrixStack {
    constructor() {
        this.transforms = [
            [
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1,
            ],
        ];
    }
    top() {
        return this.transforms[this.transforms.length - 1];
    }
    push() {
        this.transforms.push(mat4.clone(this.top()));
    }
    pop() {
        this.transforms.pop();
    }
    scale(s) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 0] *= s[0];
            m[row + 4] *= s[1];
            m[row + 8] *= s[2];
        }
    }
    scaleXYZ(x, y, z) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 0] *= x;
            m[row + 4] *= y;
            m[row + 8] *= z;
        }
    }
    scaleX(x) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 0] *= x;
        }
    }
    scaleY(y) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 4] *= y;
        }
    }
    scaleZ(z) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 8] *= z;
        }
    }
    translate(t) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 12] += m[row + 0] * t[0] + m[row + 4] * t[1] + m[row + 8] * t[2];
        }
    }
    translateXYZ(x, y, z) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 12] += m[row + 0] * x + m[row + 4] * y + m[row + 8] * z;
        }
    }
    translateX(x) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 12] += m[row + 0] * x;
        }
    }
    translateY(y) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 12] += m[row + 4] * y;
        }
    }
    translateZ(z) {
        const m = this.top();
        for (let row = 0; row < 4; ++row) {
            m[row + 12] += m[row + 8] * z;
        }
    }
    rotateX(angle) {
        const m = this.top();
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        for (let row = 0; row < 4; ++row) {
            let y = m[row + 4];
            let z = m[row + 8];
            m[row + 4] = y * c + z * s;
            m[row + 8] = z * c - y * s;
        }
    }
    rotateY(angle) {
        const m = this.top();
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        for (let row = 0; row < 4; ++row) {
            let x = m[row + 0];
            let z = m[row + 8];
            m[row + 0] = x * c - z * s;
            m[row + 8] = z * c + x * s;
        }
    }
    rotateZ(angle) {
        const m = this.top();
        const c = Math.cos(angle);
        const s = Math.sin(angle);
        for (let row = 0; row < 4; ++row) {
            let x = m[row + 0];
            let y = m[row + 4];
            m[row + 0] = x * c + y * s;
            m[row + 4] = y * c - x * s;
        }
    }
}
