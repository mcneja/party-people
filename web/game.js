"use strict";

window.onload = loadResourcesThenRun;

class Float64Grid {
    constructor(sizeX, sizeY, initialValue) {
        this.sizeX = sizeX;
        this.sizeY = sizeY;
        this.values = new Float64Array(sizeX * sizeY);
        this.fill(initialValue);
    }

    fill(value) {
        this.values.fill(value);
    }

    get(x, y) {
        return this.values[this.sizeX * y + x];
    }

    set(x, y, value) {
        this.values[this.sizeX * y + x] = value;
    }
}

class Level {
    constructor(sizeX, sizeY, initialValue) {
        this.sizeX = sizeX;
        this.sizeY = sizeY;
        this.values = new Uint8Array(sizeX * sizeY);
        this.values.fill(initialValue);
    }

    fill(value) {
        this.values.fill(value);
    }

    get(x, y) {
        return this.values[this.sizeX * y + x];
    }

    set(x, y, value) {
        this.values[this.sizeX * y + x] = value;
    }
}

function loadResourcesThenRun() {
    loadImage('font.png').then((fontImage) => { main(fontImage); });
}

function main(fontImage) {

    const canvas = document.querySelector("#canvas");
    const gl = canvas.getContext("webgl", { alpha: false, depth: false });

    if (gl == null) {
        alert("Unable to initialize WebGL. Your browser or machine may not support it.");
        return;
    }

    const renderer = createRenderer(gl, fontImage);
    const state = initState();

    canvas.requestPointerLock = canvas.requestPointerLock || canvas.mozRequestPointerLock;
    document.exitPointerLock = document.exitPointerLock || document.mozExitPointerLock;

    canvas.onclick = () => {
        if (state.paused) {
            canvas.requestPointerLock();
        } else {
            resetState(state);
            document.exitPointerLock();
        }
    };

    function requestUpdateAndRender() {
        requestAnimationFrame(now => updateAndRender(now, renderer, state));
    }

    function onLockChanged() {
        const mouseCaptured =
            document.pointerLockElement === canvas ||
            document.mozPointerLockElement === canvas;
        if (mouseCaptured) {
            document.addEventListener("mousemove", onMouseMoved, false);
            if (state.paused) {
                state.paused = false;
                state.tLast = undefined;
                state.player.velocity.x = 0;
                state.player.velocity.y = 0;
                requestUpdateAndRender();
            }
        } else {
            document.removeEventListener("mousemove", onMouseMoved, false);
            state.paused = true;
        }
    }

    function onMouseMoved(e) {
        updatePosition(state, e);
    }

    function onWindowResized() {
        requestUpdateAndRender();
    }

    document.addEventListener('pointerlockchange', onLockChanged, false);
    document.addEventListener('mozpointerlockchange', onLockChanged, false);

    window.addEventListener('resize', onWindowResized);

    requestUpdateAndRender();
}

const loadImage = src =>
    new Promise((resolve, reject) => {
        const img = new Image();
        img.onload = () => resolve(img);
        img.onerror = reject;
        img.src = src;
    });

function updatePosition(state, e) {
    if (!state.player.dead) {
        const sensitivity = 0.002;
        state.player.velocity.x += e.movementX * sensitivity;
        state.player.velocity.y -= e.movementY * sensitivity;
    }
}

function createRenderer(gl, fontImage) {
    gl.getExtension('OES_standard_derivatives');

    const renderer = {
        beginFrame: createBeginFrame(gl),
        renderField: createFieldRenderer(gl),
        renderDiscs: createDiscRenderer(gl),
        renderGlyphs: createGlyphRenderer(gl, fontImage),
        renderColoredTriangles: createColoredTriangleRenderer(gl),
    };

    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    gl.enable(gl.BLEND);
    gl.clearColor(0, 0, 0, 1);

    return renderer;
}

function initState() {
    const state = {};
    resetState(state);
    return state;
}

function resetState(state) {
    const gridSizeX = 64;
    const gridSizeY = 64;

    const player = {
        radius: 0.0125,
        position: { x: 0.5, y: 0.5 },
        velocity: { x: 0, y: 0 },
        color: { r: 0.25, g: 0.25, b: 0.25 },
        dead: false,
    };

    const enemyRadius = 0.0125;

    const obstacles = []; // createObstacles(player.position);
    const collectibles = []; // createCollectibles(obstacles);
    const costRateField = createCostRateField(gridSizeX, gridSizeY, enemyRadius, obstacles);
    const distanceFromWallsField = createDistanceFromWallsField(costRateField);
    const costRateFieldSmooth = createSmoothedCostRateField(distanceFromWallsField);
    const distanceField = createDistanceField(costRateFieldSmooth, player.position);
    const discs = []; // createEnemies(obstacles, enemyRadius, player.position);

    state.costRateField = costRateFieldSmooth;
    state.distanceField = distanceField;
    state.paused = true;
    state.gameOver = false;
    state.tLast = undefined;
    state.discs = discs;
    state.obstacles = obstacles;
    state.collectibles = collectibles;
    state.player = player;
    state.level = createLevel();
}

function createEnemies(obstacles, enemyRadius, playerPosition) {
    const enemies = [];
    const separationFromObstacle = 0.02 + enemyRadius;
    const separationFromAlly = 0.05;
    const enemyColor = { r: 0.25, g: 0.25, b: 0.25 };
    const playerDisc = { radius: 0.25, position: playerPosition };
    const angle = Math.random() * Math.PI * 2;
    const dirX = Math.cos(angle);
    const dirY = Math.sin(angle);
    for (let i = 0; i < 1000 && enemies.length < 128; ++i) {
        const enemy = {
            radius: enemyRadius,
            position: {
                x: separationFromObstacle + (1 - 2*separationFromObstacle) * Math.random(),
                y: separationFromObstacle + (1 - 2*separationFromObstacle) * Math.random(),
            },
            velocity: {
                x: 0,
                y: 0,
            },
            color: enemyColor,
            dead: false,
        };
        const dx = enemy.position.x - 0.5;
        const dy = enemy.position.y - 0.5;
        const d = dirX * dx + dirY * dy;
        if (d < 0) {
            continue;
        }
        if (!discOverlapsDiscs(enemy, obstacles, separationFromObstacle - enemyRadius) &&
            !discOverlapsDiscs(enemy, enemies, separationFromAlly) &&
            !discsOverlap(enemy, playerDisc)) {
            enemies.push(enemy);
        }
    }
    return enemies;
}

function createObstacles(playerPosition) {
    const obstacles = [];
    const radius = 0.05;
    const separation = -0.02;
    const playerDisc = { radius: 0.05, position: playerPosition };
    const color = { r: 0.25, g: 0.25, b: 0.25 };
    for (let i = 0; i < 1000 && obstacles.length < 16; ++i) {
        const obstacle = {
            radius: radius,
            position: {
                x: radius + (1 - 2*radius) * Math.random(),
                y: radius + (1 - 2*radius) * Math.random(),
            },
            color: color,
        };
        if (!discOverlapsDiscs(obstacle, obstacles, separation) &&
            !discsOverlap(obstacle, playerDisc)) {
            obstacles.push(obstacle);
        }
    }
    return obstacles;
}

function createCollectibles(obstacles) {
    const collectibles = [];
    const radius = 0.0125;
    const separationFromObstacle = 0.02;
    const separationFromCollectible = 0.05;
    const color = { r: 0.25, g: 0.25, b: 0.25 };
    for (let i = 0; i < 1000 && collectibles.length < 128; ++i) {
        const collectible = {
            radius: radius,
            position: {
                x: radius + separationFromObstacle + (1 - 2*(radius + separationFromObstacle)) * Math.random(),
                y: radius + separationFromObstacle + (1 - 2*(radius + separationFromObstacle)) * Math.random(),
            },
            color: color,
        };
        if (!discOverlapsDiscs(collectible, obstacles, separationFromObstacle) &&
            !discOverlapsDiscs(collectible, collectibles, separationFromCollectible)) {
            collectibles.push(collectible);
        }
    }
    return collectibles;
}

function discOverlapsDiscs(disc, discs, minSeparation) {
    for (const disc2 of discs) {
        const dx = disc2.position.x - disc.position.x;
        const dy = disc2.position.y - disc.position.y;
        if (dx**2 + dy**2 < (disc2.radius + disc.radius + minSeparation)**2) {
            return true;
        }
    }
    return false;
}

function discsOverlap(disc0, disc1) {
    const dx = disc1.position.x - disc0.position.x;
    const dy = disc1.position.y - disc0.position.y;
    return dx**2 + dy**2 < (disc1.radius + disc0.radius)**2;
}

function createBeginFrame(gl) {
    return () => {
        resizeCanvasToDisplaySize(gl.canvas);

        const screenX = gl.canvas.clientWidth;
        const screenY = gl.canvas.clientHeight;
    
        gl.viewport(0, 0, screenX, screenY);
        gl.clear(gl.COLOR_BUFFER_BIT);
    }
}

function createFieldRenderer(gl) {
    const vsSource = `
        attribute vec3 vPosition;
        attribute vec2 vDistance;
        attribute vec2 vSpeed;
        
        uniform mat4 uProjectionMatrix;

        varying highp float fYBlend;
        varying highp vec2 fDistance;
        varying highp vec2 fSpeed;

        void main() {
            gl_Position = uProjectionMatrix * vec4(vPosition.xy, 0, 1);
            fYBlend = vPosition.z;
            fDistance = vDistance;
            fSpeed = vSpeed;
        }
    `;

    const fsSource = `
        varying highp float fYBlend;
        varying highp vec2 fDistance;
        varying highp vec2 fSpeed;

        uniform highp float uScroll;
        uniform sampler2D uContour;

        void main() {
            highp float gamma = 2.2;
            highp float distance = mix(fDistance.x, fDistance.y, fYBlend);
            highp float speed = mix(fSpeed.x, fSpeed.y, fYBlend);
            highp vec3 speedColorLinear = vec3(1, speed, speed);
            highp float s = distance + uScroll;
            highp vec3 distanceColorLinear = pow(texture2D(uContour, vec2(s, 0)).rgb, vec3(gamma));
            highp vec3 colorLinear = speedColorLinear * distanceColorLinear;
            gl_FragColor.rgb = pow(colorLinear, vec3(1.0/gamma));
            gl_FragColor.a = 1.0;
        }
    `;

    const projectionMatrix = new Float32Array(16);
    projectionMatrix.fill(0);
    projectionMatrix[0] = 1;
    projectionMatrix[5] = 1;
    projectionMatrix[10] = 1;
    projectionMatrix[12] = -1;
    projectionMatrix[13] = -1;
    projectionMatrix[15] = 1;

    const program = initShaderProgram(gl, vsSource, fsSource);

    const vertexAttributeLoc = {
        position: gl.getAttribLocation(program, 'vPosition'),
        distance: gl.getAttribLocation(program, 'vDistance'),
        speed: gl.getAttribLocation(program, 'vSpeed'),
    };

    const uniformLoc = {
        projectionMatrix: gl.getUniformLocation(program, 'uProjectionMatrix'),
        uScroll: gl.getUniformLocation(program, 'uScroll'),
        uContour: gl.getUniformLocation(program, 'uContour'),
    };

    const vertexBuffer = gl.createBuffer();

    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    const stride = 28; // seven 4-byte floats
    gl.vertexAttribPointer(vertexAttributeLoc.position, 3, gl.FLOAT, false, stride, 0);
    gl.vertexAttribPointer(vertexAttributeLoc.distance, 2, gl.FLOAT, false, stride, 12);
    gl.vertexAttribPointer(vertexAttributeLoc.speed, 2, gl.FLOAT, false, stride, 20);

    const contourTexture = createStripeTexture(gl);

    return (costRateField, distanceField, uScroll) => {
        gl.useProgram(program);

        gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);

        const vertexInfo = createVertexInfo(costRateField, distanceField);
        gl.bufferData(gl.ARRAY_BUFFER, vertexInfo, gl.DYNAMIC_DRAW);

        gl.enableVertexAttribArray(vertexAttributeLoc.position);
        gl.enableVertexAttribArray(vertexAttributeLoc.distance);
        gl.enableVertexAttribArray(vertexAttributeLoc.speed);
        const stride = 28; // seven 4-byte floats
        gl.vertexAttribPointer(vertexAttributeLoc.position, 3, gl.FLOAT, false, stride, 0);
        gl.vertexAttribPointer(vertexAttributeLoc.distance, 2, gl.FLOAT, false, stride, 12);
        gl.vertexAttribPointer(vertexAttributeLoc.speed, 2, gl.FLOAT, false, stride, 20);
    
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, contourTexture);
    
        gl.uniform1i(uniformLoc.uContour, 0);

        const gridSizeX = costRateField.sizeX;
        const gridSizeY = costRateField.sizeY;

        projectionMatrix[0] = 2 / (gridSizeX - 1);
        projectionMatrix[5] = 2 / (gridSizeY - 1);
        gl.uniformMatrix4fv(uniformLoc.projectionMatrix, false, projectionMatrix);
    
        gl.uniform1f(uniformLoc.uScroll, uScroll);
    
        gl.drawArrays(gl.TRIANGLES, 0, (gridSizeX - 1) * (gridSizeY - 1) * 6);
    
        gl.disableVertexAttribArray(vertexAttributeLoc.speed);
        gl.disableVertexAttribArray(vertexAttributeLoc.distance);
        gl.disableVertexAttribArray(vertexAttributeLoc.position);
    };
}

function createDiscRenderer(gl) {
    const vsSource = `
        attribute vec2 vPosition;
        
        uniform mat4 uProjectionMatrix;

        varying highp vec2 fPosition;

        void main() {
            fPosition = vPosition;
            gl_Position = uProjectionMatrix * vec4(vPosition.xy, 0, 1);
        }
    `;

    const fsSource = `
        #extension GL_OES_standard_derivatives : enable

        varying highp vec2 fPosition;

        uniform highp vec3 uColor;

        void main() {
            highp float r = length(fPosition);
            highp float aaf = fwidth(r);
            highp float opacity = 1.0 - smoothstep(1.0 - aaf, 1.0, r);
            gl_FragColor = vec4(uColor, opacity);
        }
    `;

    const projectionMatrix = new Float32Array(16);
    projectionMatrix.fill(0);
    projectionMatrix[10] = 1;
    projectionMatrix[15] = 1;

    const program = initShaderProgram(gl, vsSource, fsSource);

    const vertexPositionLoc = gl.getAttribLocation(program, 'vPosition');
    const projectionMatrixLoc = gl.getUniformLocation(program, 'uProjectionMatrix');
    const colorLoc = gl.getUniformLocation(program, 'uColor');
    const vertexBuffer = createDiscVertexBuffer(gl);

    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    const stride = 8; // two 4-byte floats
    gl.vertexAttribPointer(vertexPositionLoc, 2, gl.FLOAT, false, stride, 0);

    return discs => {
        gl.useProgram(program);

        gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
        gl.enableVertexAttribArray(vertexPositionLoc);
        gl.vertexAttribPointer(vertexPositionLoc, 2, gl.FLOAT, false, 0, 0);

        for (const disc of discs) {
            gl.uniform3f(colorLoc, disc.color.r, disc.color.g, disc.color.b);

            projectionMatrix[0] = 2 * disc.radius;
            projectionMatrix[5] = 2 * disc.radius;
            projectionMatrix[12] = 2 * disc.position.x - 1;
            projectionMatrix[13] = 2 * disc.position.y - 1;
            gl.uniformMatrix4fv(projectionMatrixLoc, false, projectionMatrix);

            gl.drawArrays(gl.TRIANGLES, 0, 6);
        }
    
        gl.disableVertexAttribArray(vertexPositionLoc);
    };
}

function createDiscVertexBuffer(gl) {
    const v = new Float32Array(6 * 2);
    let i = 0;

    function makeVert(x, y) {
        v[i++] = x;
        v[i++] = y;
    }

    makeVert(-1, -1);
    makeVert( 1, -1);
    makeVert( 1,  1);
    makeVert( 1,  1);
    makeVert(-1,  1);
    makeVert(-1, -1);

    const vertexBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, v, gl.STATIC_DRAW);

    return vertexBuffer;
}

function createVertexInfo(costRateField, distanceField) {
    const sizeX = costRateField.sizeX;
    const sizeY = costRateField.sizeY;
    const v = new Float32Array(7 * 6 * (sizeX - 1) * (sizeY - 1));
    let i = 0;

    function distance(x, y) {
        return distanceField.get(x, y);
    }

    function speed(x, y) {
        return 1 / costRateField.get(x, y);
    }

    function makeVert(x, y, s, d0, d1, c0, c1) {
        v[i++] = x;
        v[i++] = y;
        v[i++] = s;
        v[i++] = d0;
        v[i++] = d1;
        v[i++] = c0;
        v[i++] = c1;
    }

    for (let y = 0; y < sizeY - 1; ++y) {
        for (let x = 0; x < sizeX - 1; ++x) {
            const dist00 = distance(x, y);
            const dist10 = distance(x+1, y);
            const dist01 = distance(x, y+1);
            const dist11 = distance(x+1, y+1);
            const speed00 = speed(x, y);
            const speed10 = speed(x+1, y);
            const speed01 = speed(x, y+1);
            const speed11 = speed(x+1, y+1);
            makeVert(x, y, 0, dist00, dist01, speed00, speed01);
            makeVert(x+1, y, 0, dist10, dist11, speed10, speed11);
            makeVert(x, y+1, 1, dist00, dist01, speed00, speed01);
            makeVert(x, y+1, 1, dist00, dist01, speed00, speed01);
            makeVert(x+1, y, 0, dist10, dist11, speed10, speed11);
            makeVert(x+1, y+1, 1, dist10, dist11, speed10, speed11);
        }
    }

    return v;
}

function createGlyphRenderer(gl, fontImage) {
    const vsSource = `
        attribute vec4 vPositionTexcoord;
        attribute vec4 vColor;

        uniform mat4 uProjectionMatrix;

        varying highp vec2 fTexcoord;
        varying highp vec4 fColor;

        void main() {
            fTexcoord = vPositionTexcoord.zw;
            fColor = vColor;
            gl_Position = uProjectionMatrix * vec4(vPositionTexcoord.xy, 0, 1);
        }
    `;

    const fsSource = `
        varying highp vec2 fTexcoord;
        varying highp vec4 fColor;

        uniform sampler2D uOpacity;

        void main() {
            gl_FragColor = fColor * vec4(1, 1, 1, texture2D(uOpacity, fTexcoord));
        }
    `;

    const fontTexture = createTextureFromImage(gl, fontImage);

    const projectionMatrixData = new Float32Array(16);
    projectionMatrixData.fill(0);
    projectionMatrixData[0] = 1;
    projectionMatrixData[5] = 1;
    projectionMatrixData[10] = 1;
    projectionMatrixData[15] = 1;

    const program = initShaderProgram(gl, vsSource, fsSource);

    const vPositionTexcoordLoc = gl.getAttribLocation(program, 'vPositionTexcoord');
    const vColorLoc = gl.getAttribLocation(program, 'vColor');

    const uProjectionMatrixLoc = gl.getUniformLocation(program, 'uProjectionMatrix');
    const uOpacityLoc = gl.getUniformLocation(program, 'uOpacity');

    const maxQuads = 4096;
    const numVertices = 4 * maxQuads;
    const bytesPerVertex = 4 * Float32Array.BYTES_PER_ELEMENT + Uint32Array.BYTES_PER_ELEMENT;
    const wordsPerQuad = bytesPerVertex;

    const indexBuffer = createGlyphIndexBuffer(gl, maxQuads);

    const vertexBuffer = gl.createBuffer();

    const vertexData = new ArrayBuffer(numVertices * bytesPerVertex);
    const vertexDataAsFloat32 = new Float32Array(vertexData);
    const vertexDataAsUint32 = new Uint32Array(vertexData);

    let numQuads = 0;

    function addQuad(x0, y0, x1, y1, s0, t0, s1, t1, color) {
        if (numQuads >= maxQuads) {
            flushQuads();
        }

        let i = numQuads * wordsPerQuad;

        vertexDataAsFloat32[i+0] = x0;
        vertexDataAsFloat32[i+1] = y0;
        vertexDataAsFloat32[i+2] = s0;
        vertexDataAsFloat32[i+3] = t0;
        vertexDataAsUint32[i+4] = color;

        vertexDataAsFloat32[i+5] = x1;
        vertexDataAsFloat32[i+6] = y0;
        vertexDataAsFloat32[i+7] = s1;
        vertexDataAsFloat32[i+8] = t0;
        vertexDataAsUint32[i+9] = color;

        vertexDataAsFloat32[i+10] = x0;
        vertexDataAsFloat32[i+11] = y1;
        vertexDataAsFloat32[i+12] = s0;
        vertexDataAsFloat32[i+13] = t1;
        vertexDataAsUint32[i+14] = color;

        vertexDataAsFloat32[i+15] = x1;
        vertexDataAsFloat32[i+16] = y1;
        vertexDataAsFloat32[i+17] = s1;
        vertexDataAsFloat32[i+18] = t1;
        vertexDataAsUint32[i+19] = color;

        ++numQuads;
    }

    function flushQuads() {
        if (numQuads <= 0) {
            return;
        }

        gl.useProgram(program);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, fontTexture);
        gl.uniform1i(uOpacityLoc, 0);

        gl.uniformMatrix4fv(uProjectionMatrixLoc, false, projectionMatrixData);

        gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, vertexData, gl.DYNAMIC_DRAW);
        gl.vertexAttribPointer(vPositionTexcoordLoc, 4, gl.FLOAT, false, bytesPerVertex, 0);
        gl.vertexAttribPointer(vColorLoc, 4, gl.UNSIGNED_BYTE, true, bytesPerVertex, 16);
        gl.enableVertexAttribArray(vPositionTexcoordLoc);
        gl.enableVertexAttribArray(vColorLoc);

        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer);

        gl.drawElements(gl.TRIANGLES, 6 * numQuads, gl.UNSIGNED_SHORT, 0);

        numQuads = 0;
    }

    return {
        add: addQuad,
        flush: flushQuads,
    };
}

function createGlyphIndexBuffer(gl, maxQuads) {
    const indices = new Uint16Array(maxQuads * 6);

    for (let i = 0; i < maxQuads; ++i) {
        let j = 6*i;
        let k = 4*i;
        indices[j+0] = k+0;
        indices[j+1] = k+1;
        indices[j+2] = k+2;
        indices[j+3] = k+2;
        indices[j+4] = k+1;
        indices[j+5] = k+3;
    }

    const indexBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indices, gl.STATIC_DRAW);

    return indexBuffer;
}

function createTextureFromImage(gl, image) {
    const level = 0;
    const internalFormat = gl.RGBA;
    const srcFormat = gl.RGBA;
    const srcType = gl.UNSIGNED_BYTE;
    const texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texImage2D(gl.TEXTURE_2D, level, internalFormat, srcFormat, srcType, image);
    gl.generateMipmap(gl.TEXTURE_2D);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    return texture;
}

function updateAndRender(now, renderer, state) {
    const t = now / 1000;
    const dt = (state.paused || state.tLast === undefined) ? 0 : Math.min(1/30, t - state.tLast);
    state.tLast = t;

    if (dt > 0) {
        updateState(state, dt);
    }

    drawScreen(renderer, state);

    if (!state.paused) {
        requestAnimationFrame(now => updateAndRender(now, renderer, state));
    }
}

function createColoredTriangleRenderer(gl) {
    const vsSource = `
        attribute vec2 vPosition;
        attribute vec4 vColor;

        uniform mat4 uProjectionMatrix;

        varying highp vec4 fColor;

        void main() {
            fColor = vColor;
            gl_Position = uProjectionMatrix * vec4(vPosition.xy, 0, 1);
        }
    `;

    const fsSource = `
        varying highp vec4 fColor;
        void main() {
            gl_FragColor = fColor;
        }
    `;

    const projectionMatrix = new Float32Array(16);
    projectionMatrix.fill(0);
    projectionMatrix[10] = 1;
    projectionMatrix[12] = -1;
    projectionMatrix[13] = -1;
    projectionMatrix[15] = 1;

    const program = initShaderProgram(gl, vsSource, fsSource);

    const vertexPositionLoc = gl.getAttribLocation(program, 'vPosition');
    const vertexColorLoc = gl.getAttribLocation(program, 'vColor');
    const projectionMatrixLoc = gl.getUniformLocation(program, 'uProjectionMatrix');

    const vertexBuffer = gl.createBuffer();

    const bytesPerVertex = 12; // two 4-byte floats and one 32-bit color
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.vertexAttribPointer(vertexPositionLoc, 2, gl.FLOAT, false, bytesPerVertex, 0);
    gl.vertexAttribPointer(vertexColorLoc, 4, gl.UNSIGNED_BYTE, true, bytesPerVertex, 8);

    return (worldSizeX, worldSizeY, vertexData) => {
        const numVerts = Math.floor(vertexData.byteLength / bytesPerVertex);

        gl.useProgram(program);

        projectionMatrix[0] = 2 / worldSizeX;
        projectionMatrix[5] = 2 / worldSizeY;
        gl.uniformMatrix4fv(projectionMatrixLoc, false, projectionMatrix);

        gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, vertexData, gl.DYNAMIC_DRAW);
        gl.enableVertexAttribArray(vertexPositionLoc);
        gl.enableVertexAttribArray(vertexColorLoc);
        gl.vertexAttribPointer(vertexPositionLoc, 2, gl.FLOAT, false, bytesPerVertex, 0);
        gl.vertexAttribPointer(vertexColorLoc, 4, gl.UNSIGNED_BYTE, true, bytesPerVertex, 8);
    
        gl.drawArrays(gl.TRIANGLES, 0, numVerts);

        gl.disableVertexAttribArray(vertexColorLoc);
        gl.disableVertexAttribArray(vertexPositionLoc);
    };
}

function updateState(state, dt) {

    if (state.player.dead) {
        const r = Math.exp(-dt);
        state.player.velocity.x *= r;
        state.player.velocity.y *= r;
    }

    state.player.position.x += state.player.velocity.x * dt;
    state.player.position.y += state.player.velocity.y * dt;

    for (const obstacle of state.obstacles) {
        fixupPositionAndVelocityAgainstDisc(state.player, obstacle);
    }

    fixupPositionAndVelocityAgainstBoundary(state.player);

    if (!state.gameOver) {
        state.collectibles = state.collectibles.filter(collectible => !discsOverlap(state.player, collectible));
    }
    if (state.collectibles.length <= 0) {
        state.gameOver = true;
    }

    updateDistanceField(state.costRateField, state.distanceField, state.player.position);

    const discSpeed = 0.2 - 0.15 * Math.min(state.collectibles.length, 80) / 80;

    let enemyDied = false;
    for (const disc of state.discs) {
        updateEnemy(state.distanceField, dt, discSpeed, state.player, disc);
        if (disc.dead) {
            enemyDied = true;
            if (!state.gameOver) {
                state.player.dead = true;
                state.player.color.r *= 0.5;
                state.player.color.g *= 0.5;
                state.player.color.b *= 0.5;
                state.gameOver = true;
            }
        }
    }

    if (enemyDied) {
        state.discs = state.discs.filter(disc => !disc.dead);
    }

    for (let k = 0; k < 3; ++k) {
        for (let i = 0; i < state.discs.length; ++i) {
            for (let j = i + 1; j < state.discs.length; ++j) {
                fixupDiscPairs(state.discs[i], state.discs[j]);
            }

            for (const obstacle of state.obstacles) {
                fixupPositionAndVelocityAgainstDisc(state.discs[i], obstacle);
            }

            fixupPositionAndVelocityAgainstBoundary(state.discs[i]);
        }
    }
}

function updateEnemy(distanceField, dt, discSpeed, player, disc) {
    if (discsOverlap(disc, player)) {
        disc.dead = true;
        return;
    }

    const gradient = estimateGradient(distanceField, disc.position.x, disc.position.y);

    const gradientLen = Math.sqrt(gradient.x**2 + gradient.y**2);

    const scale = discSpeed / Math.max(1e-8, gradientLen);

    const vxPrev = disc.velocity.x;
    const vyPrev = disc.velocity.y;

    const vxTarget = gradient.x * -scale;
    const vyTarget = gradient.y * -scale;

    let dvx = vxTarget - disc.velocity.x;
    let dvy = vyTarget - disc.velocity.y;
    const dv = Math.sqrt(dvx**2 + dvy**2);

    const maxDv = 25 * (discSpeed**2) * dt;

    if (dv > maxDv) {
        dvx *= maxDv / dv;
        dvy *= maxDv / dv;
    }

    disc.velocity.x += dvx;
    disc.velocity.y += dvy;

    disc.position.x += (vxPrev + disc.velocity.x) * dt / 2;
    disc.position.y += (vyPrev + disc.velocity.y) * dt / 2;
}

function fixupDiscPairs(disc0, disc1) {
    const dx = disc1.position.x - disc0.position.x;
    const dy = disc1.position.y - disc0.position.y;
    const d = Math.sqrt(dx**2 + dy**2);
    const dist = d - (disc0.radius + disc1.radius);

    if (dist < 0) {
        disc0.position.x += dx * dist / (2 * d);
        disc0.position.y += dy * dist / (2 * d);
        disc1.position.x -= dx * dist / (2 * d);
        disc1.position.y -= dy * dist / (2 * d);

        const vx = disc1.velocity.x - disc0.velocity.x;
        const vy = disc1.velocity.y - disc0.velocity.y;
        const vn = vx * dx + vy * dy;
        if (vn < 0) {
            disc0.velocity.x += vn * dx / (2 * d**2);
            disc0.velocity.y += vn * dy / (2 * d**2);
            disc1.velocity.x -= vn * dx / (2 * d**2);
            disc1.velocity.y -= vn * dy / (2 * d**2);
        }
    }
}

function fixupPositionAndVelocityAgainstBoundary(disc) {
    if (disc.position.x < disc.radius) {
        disc.position.x = disc.radius;
        disc.velocity.x = 0;
    } else if (disc.position.x > 1 - disc.radius) {
        disc.position.x = 1 - disc.radius;
        disc.velocity.x = 0;
    }
    if (disc.position.y < disc.radius) {
        disc.position.y = disc.radius;
        disc.velocity.y = 0;
    } else if (disc.position.y > 1 - disc.radius) {
        disc.position.y = 1 - disc.radius;
        disc.velocity.y = 0;
    }
}

function fixupPositionAndVelocityAgainstDisc(disc, obstacle) {
    const dx = disc.position.x - obstacle.position.x;
    const dy = disc.position.y - obstacle.position.y;
    const d = Math.sqrt(dx**2 + dy**2);
    const dist = d - (disc.radius + obstacle.radius);

    if (dist < 0) {
        disc.position.x -= dx * dist / d;
        disc.position.y -= dy * dist / d;

        const vn = disc.velocity.x * dx + disc.velocity.y * dy;
        if (vn < 0) {
            disc.velocity.x -= vn * dx / d**2;
            disc.velocity.y -= vn * dy / d**2;
        }
    }
}

function drawScreen(renderer, state) {
    renderer.beginFrame();
    renderer.renderColoredTriangles(state.level.worldSizeX, state.level.worldSizeY, state.level.vertexData);
//    renderer.renderField(state.costRateField, state.distanceField, 0);
//    renderer.renderDiscs(state.obstacles);
//    renderer.renderDiscs(state.collectibles);
//    renderer.renderDiscs(state.discs);
//    renderer.renderDiscs([state.player]);

    const r = 0.0125 * 1.2;
    const rx = r;
    const ry = r * 2;
    const tx = 0.0625;
    const ty = 0.0625;

    for (const c of state.collectibles) {
        const x = c.position.x * 2 - 1;
        const y = c.position.y * 2 - 1;
        renderer.renderGlyphs.add(x - rx, y - ry, x + rx, y + ry, 15*tx, ty, 16*tx, 0, 0xff00ffff);
    }

    for (const d of state.discs) {
        const x = d.position.x * 2 - 1;
        const y = d.position.y * 2 - 0.995;
        renderer.renderGlyphs.add(x - rx, y - ry, x + rx, y + ry, 7*tx, 7*ty, 8*tx, 6*ty, 0xff00a000);
    }

//    renderer.renderGlyphs.add(-0.25, -0.5, 0.25, 0.5, 0, 1, 1, 0, 0xffffffff);
    const x = state.player.position.x * 2 - 1;
    const y = state.player.position.y * 2 - 1;
    renderer.renderGlyphs.add(x - rx, y - ry, x + rx, y + ry, 1*tx, ty, 2*tx, 0, 0xff00ffff);
    renderer.renderGlyphs.flush();
}

function resizeCanvasToDisplaySize(canvas) {
    const rect = canvas.parentNode.getBoundingClientRect();
    if (canvas.width !== rect.width || canvas.height !== rect.height) {
        canvas.width = rect.width;
        canvas.height = rect.height;
    }
}

function initShaderProgram(gl, vsSource, fsSource) {
    const vertexShader = loadShader(gl, gl.VERTEX_SHADER, vsSource);
    const fragmentShader = loadShader(gl, gl.FRAGMENT_SHADER, fsSource);

    const program = gl.createProgram();
    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);

    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        alert('Unable to initialize the shader program: ' + gl.getProgramInfoLog(program));
        return null;
    }

    return program;
}

function loadShader(gl, type, source) {
    const shader = gl.createShader(type);

    gl.shaderSource(shader, source);
    gl.compileShader(shader);

    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
        alert('An error occurred compiling the shaders: ' + gl.getShaderInfoLog(shader));
        gl.deleteShader(shader);
        return null;
    }

    return shader;
}

function createStripeTexture(gl) {
    const stripeImageWidth = 64;
    const stripeImage = new Uint8Array(stripeImageWidth);
    for (let j = 0; j < stripeImageWidth; ++j) {
        stripeImage[j] = 224 + (stripeImageWidth - j) / 4;
    }

    const texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);

    const level = 0;
    const internalFormat = gl.LUMINANCE;
    const srcFormat = gl.LUMINANCE;
    const srcType = gl.UNSIGNED_BYTE;
    gl.texImage2D(gl.TEXTURE_2D, level, internalFormat, stripeImageWidth, 1, 0, srcFormat, srcType, stripeImage);
    gl.generateMipmap(gl.TEXTURE_2D);

    return texture;
}

function priorityQueuePop(q) {
    const x = q[0];
    q[0] = q[q.length - 1]; // q.at(-1);
    q.pop();
    let i = 0;
    const c = q.length;
    while (true) {
        let iChild = i;
        const iChild0 = 2*i + 1;
        if (iChild0 < c && q[iChild0].priority < q[iChild].priority) {
            iChild = iChild0;
        }
        const iChild1 = iChild0 + 1;
        if (iChild1 < c && q[iChild1].priority < q[iChild].priority) {
            iChild = iChild1;
        }
        if (iChild == i) {
            break;
        }
        [q[i], q[iChild]] = [q[iChild], q[i]];
        i = iChild;
    }
    return x;
}

function priorityQueuePush(q, x) {
    q.push(x);
    let i = q.length - 1;
    while (i > 0) {
        const iParent = Math.floor((i - 1) / 2);
        if (q[i].priority >= q[iParent].priority) {
            break;
        }
        [q[i], q[iParent]] = [q[iParent], q[i]];
        i = iParent;
    }
}

function createCostRateField(sizeX, sizeY, enemyRadius, obstacles) {
    const costRate = new Float64Grid(sizeX, sizeY, 1);

    function dist(x, y) {
        let minDist = Infinity;
        x /= sizeX - 1;
        y /= sizeY - 1;
        for (const obstacle of obstacles) {
            const dx = x - obstacle.position.x;
            const dy = y - obstacle.position.y;
            const dist = Math.sqrt(dx**2 + dy**2) - (obstacle.radius + enemyRadius);
            minDist = Math.min(minDist, dist);
        }
        minDist = Math.max(0, minDist);
        return minDist * (sizeX - 1);
    }

    for (let y = 0; y < sizeY; ++y) {
        for (let x = 0; x < sizeX; ++x) {
            const d = dist(x, y);
            costRate.set(x, y, (d <= 0) ? 1000 : 1);
        }
    }

    return costRate;
}

function createSmoothedCostRateField(distanceFromWallsField) {
    const sizeX = distanceFromWallsField.sizeX;
    const sizeY = distanceFromWallsField.sizeY;

    const costRateFieldSmooth = new Float64Grid(sizeX, sizeY, 1);

    for (let y = 0; y < sizeY; ++y) {
        for (let x = 0; x < sizeX; ++x) {
            const distance = distanceFromWallsField.get(x, y);

            const costRate = 1 + Math.min(1e6, 10 / distance**2);

            costRateFieldSmooth.set(x, y, costRate);
        }
    }
    return costRateFieldSmooth;
}

function createDistanceFromWallsField(costRateField) {
    const sizeX = costRateField.sizeX;
    const sizeY = costRateField.sizeY;
    const toVisit = [];
    for (let y = 0; y < sizeY; ++y) {
        for (let x = 0; x < sizeX; ++x) {
            if (costRateField.get(x, y) > 1.0) {
                toVisit.push({priority: 0, x: x, y: y});
            }
        }
    }

    const distanceField = new Float64Grid(sizeX, sizeY, Infinity);
    fastMarchFill(distanceField, toVisit, (x, y) => estimatedDistance(distanceField, x, y));

    const scale = 128 / sizeX;
    distanceField.values.forEach((x, i, d) => d[i] *= scale);

    return distanceField;
}

function createDistanceField(costRateField, goal) {
    const sizeX = costRateField.sizeX;
    const sizeY = costRateField.sizeY;
    const distanceField = new Float64Grid(costRateField.sizeX, costRateField.sizeY, Infinity);
    updateDistanceField(costRateField, distanceField, goal);
    return distanceField;
}

function updateDistanceField(costRateField, distanceField, goal) {
    const sizeX = costRateField.sizeX;
    const sizeY = costRateField.sizeY;
    const goalX = goal.x * (sizeX - 1);
    const goalY = goal.y * (sizeY - 1);
    const r = 2;
    const xMin = Math.min(sizeX - 1, Math.max(0, Math.floor(goalX + 0.5) - r));
    const yMin = Math.min(sizeY - 1, Math.max(0, Math.floor(goalY + 0.5) - r));
    const xMax = Math.min(sizeX, xMin + 2*r+1);
    const yMax = Math.min(sizeY, yMin + 2*r+1);

    distanceField.fill(Infinity);

    for (let y = yMin; y < yMax; ++y) {
        for (let x = xMin; x < xMax; ++x) {
            const dx = x - goalX;
            const dy = y - goalY;
            const dist = Math.sqrt(dx**2 + dy**2);
            const costRate = costRateField.get(x, y);
            const cost = dist * costRate;
            distanceField.set(x, y, cost);
        }
    }

    let toVisit = [];
    for (let y = Math.max(0, yMin - 1); y < Math.min(sizeY, yMax + 1); ++y) {
        for (let x = Math.max(0, xMin - 1); x < Math.min(sizeX, xMax + 1); ++x) {
            if (distanceField.get(x, y) === Infinity) {
                const cost = estimatedDistanceWithSpeed(distanceField, costRateField, x, y);
                priorityQueuePush(toVisit, {priority: cost, x: x, y: y});
            }
        }
    }

    fastMarchFill(distanceField, toVisit, (x, y) => estimatedDistanceWithSpeed(distanceField, costRateField, x, y));

    const scale = 32 / sizeX;
    distanceField.values.forEach((x, i, d) => d[i] *= scale);
}

function fastMarchFill(field, toVisit, estimatedDistance) {
    while (toVisit.length > 0) {
        const {priority, x, y} = priorityQueuePop(toVisit);

        if (field.get(x, y) <= priority) {
            continue;
        }

        field.set(x, y, priority);

        if (x < field.sizeX - 1) {
            const d = estimatedDistance(x + 1, y);
            if (d < field.get(x+1, y)) {
                priorityQueuePush(toVisit, {priority: d, x: x+1, y: y});
            }
        }

        if (x > 0) {
            const d = estimatedDistance(x - 1, y);
            if (d < field.get(x-1, y)) {
                priorityQueuePush(toVisit, {priority: d, x: x-1, y: y});
            }
        }

        if (y < field.sizeY - 1) {
            const d = estimatedDistance(x, y + 1);
            if (d < field.get(x, y+1)) {
                priorityQueuePush(toVisit, {priority: d, x: x, y: y+1});
            }
        }

        if (y > 0) {
            const d = estimatedDistance(x, y - 1);
            if (d < field.get(x, y-1)) {
                priorityQueuePush(toVisit, {priority: d, x: x, y: y-1});
            }
        }
    }
}

function estimatedDistance(field, x, y) {
    const dXNeg = (x > 0) ? field.get(x-1, y) : Infinity;
    const dXPos = (x < field.sizeX - 1) ? field.get(x+1, y) : Infinity;
    const dYNeg = (y > 0) ? field.get(x, y-1) : Infinity;
    const dYPos = (y < field.sizeY - 1) ? field.get(x, y+1) : Infinity;

    const dXMin = Math.min(dXNeg, dXPos);
    const dYMin = Math.min(dYNeg, dYPos);

    const timeHorizontal = 1.0;

    const d = (Math.abs(dXMin - dYMin) <= timeHorizontal) ?
        ((dXMin + dYMin) + Math.sqrt((dXMin + dYMin)**2 - 2 * (dXMin**2 + dYMin**2 - timeHorizontal**2))) / 2:
        Math.min(dXMin, dYMin) + timeHorizontal;

    return d;
}

function estimatedDistanceWithSpeed(field, speed, x, y) {
    const dXNeg = (x > 0) ? field.get(x-1, y) : Infinity;
    const dXPos = (x < field.sizeX - 1) ? field.get(x+1, y) : Infinity;
    const dYNeg = (y > 0) ? field.get(x, y-1) : Infinity;
    const dYPos = (y < field.sizeY - 1) ? field.get(x, y+1) : Infinity;

    const dXMin = Math.min(dXNeg, dXPos);
    const dYMin = Math.min(dYNeg, dYPos);

    const timeHorizontal = speed.get(x, y);

    const d = (Math.abs(dXMin - dYMin) <= timeHorizontal) ?
        ((dXMin + dYMin) + Math.sqrt((dXMin + dYMin)**2 - 2 * (dXMin**2 + dYMin**2 - timeHorizontal**2))) / 2:
        Math.min(dXMin, dYMin) + timeHorizontal;

    return d;
}

function safeDistance(distanceField, x, y) {
    if (x < 0 || y < 0 || x >= distanceField.sizeX || y >= distanceField.sizeY)
        return 1e4;

    return distanceField.get(x, y);
}

function estimateGradient(distanceField, x, y) {

    x *= (distanceField.sizeX - 1);
    y *= (distanceField.sizeY - 1);

    const uX = x - Math.floor(x);
    const uY = y - Math.floor(y);

    const gridX = Math.floor(x);
    const gridY = Math.floor(y);

    const d00 = safeDistance(distanceField, gridX, gridY);
    const d10 = safeDistance(distanceField, gridX + 1, gridY);
    const d01 = safeDistance(distanceField, gridX, gridY + 1);
    const d11 = safeDistance(distanceField, gridX + 1, gridY + 1);

    const gradX = (d10 - d00) * (1 - uY) + (d11 - d01) * uY;
    const gradY = (d01 - d00) * (1 - uX) + (d11 - d10) * uX;

    return { x: gradX, y: gradY };
}

function estimateDistance(distanceField, x, y) {

    x *= (distanceField.sizeX - 1);
    y *= (distanceField.sizeY - 1);

    const uX = x - Math.floor(x);
    const uY = y - Math.floor(y);

    const gridX = Math.floor(x);
    const gridY = Math.floor(y);

    const d00 = safeDistance(distanceField, gridX, gridY);
    const d10 = safeDistance(distanceField, gridX + 1, gridY);
    const d01 = safeDistance(distanceField, gridX, gridY + 1);
    const d11 = safeDistance(distanceField, gridX + 1, gridY + 1);

    const d0 = d00 + (d10 - d00) * uX;
    const d1 = d01 + (d11 - d01) * uX;
    const d = d0 + (d1 - d0) * uY;

    return d;
}

function randomInRange(n) {
    return Math.floor(Math.random() * n);
}

// The output level data structure consists of dimensions and
// a byte array with the tile for each square. The starting
// position is also returned.

function coinFlips(total) {
    let count = 0;
    while (total > 0) {
        if (Math.random() < 0.5)
            ++count;
        --total;
    }
    return count;
}

const ttSolid = 0;
const ttRoom = 1;
const ttHall = 2;
const ttWall = 3;

function createLevel() {
    const mapSizeX = 96;
    const mapSizeY = 72;

    const numCellsX = 4;
    const numCellsY = 3;
    const corridorWidth = 5;

    const squaresPerBlockX = Math.floor((mapSizeX + corridorWidth) / numCellsX);
    const squaresPerBlockY = Math.floor((mapSizeY + corridorWidth) / numCellsY);

    const minRoomSize = corridorWidth + 2;

    // Create some rooms.

    const rooms = [];

    for (let roomY = 0; roomY < numCellsY; ++roomY) {
        for (let roomX = 0; roomX < numCellsX; ++roomX) {
            const cellMinX = roomX * squaresPerBlockX + 1;
            const cellMinY = roomY * squaresPerBlockY + 1;
            const cellMaxX = (roomX + 1) * squaresPerBlockX - (1 + corridorWidth);
            const cellMaxY = (roomY + 1) * squaresPerBlockY - (1 + corridorWidth);
            const maxRoomSizeX = cellMaxX - cellMinX;
            const maxRoomSizeY = cellMaxY - cellMinY;
            const halfRoomSizeRangeX = Math.floor((maxRoomSizeX + 1 - minRoomSize) / 2);
            const halfRoomSizeRangeY = Math.floor((maxRoomSizeY + 1 - minRoomSize) / 2);

            const roomSizeX = 2 * coinFlips(halfRoomSizeRangeX) + minRoomSize;
            const roomSizeY = 2 * coinFlips(halfRoomSizeRangeY) + minRoomSize;

            const roomMinX = randomInRange(1 + maxRoomSizeX - roomSizeX) + cellMinX;
            const roomMinY = randomInRange(1 + maxRoomSizeY - roomSizeY) + cellMinY;

            rooms.push({
                minX: roomMinX,
                minY: roomMinY,
                sizeX: roomSizeX,
                sizeY: roomSizeY,
            });
        }
    }

    // Generate the graph of connections between rooms

    const potentialEdges = [];
    for (let roomY = 0; roomY < numCellsY; ++roomY) {
        for (let roomX = 1; roomX < numCellsX; ++roomX) {
            const room1 = roomY * numCellsX + roomX;
            const room0 = room1 - 1;
            potentialEdges.push([room0, room1]);
        }
    }

    for (let roomY = 1; roomY < numCellsY; ++roomY) {
        for (let roomX = 0; roomX < numCellsX; ++roomX) {
            const room1 = roomY * numCellsX + roomX;
            const room0 = room1 - numCellsX;
            potentialEdges.push([room0, room1]);
        }
    }

    shuffleArray(potentialEdges);

    const numRooms = numCellsX * numCellsY;
    const roomGroup = [];
    for (let i = 0; i < numRooms; ++i) {
        roomGroup.push(i);
    }

    const edges = [];
    for (const edge of potentialEdges) {
        const group0 = roomGroup[edge[0]];
        const group1 = roomGroup[edge[1]];
        const span = group0 != group1;
        if (span || Math.random() < 0.5) {
            edges.push({edge: edge, span: span});
            for (let i = 0; i < numRooms; ++i) {
                if (roomGroup[i] === group1) {
                    roomGroup[i] = group0;
                }
            }
        }
    }

    // Plot rooms into a grid

    const level = new Level(mapSizeX, mapSizeY, ttSolid);

    for (const room of rooms) {
        for (let y = 0; y < room.sizeY; ++y) {
            for (let x = 0; x < room.sizeX; ++x) {
                level.set(x + room.minX, y + room.minY, ttRoom);
            }
        }

        for (let x = 0; x < room.sizeX; ++x) {
            level.set(x + room.minX, room.minY - 1, ttWall);
            level.set(x + room.minX, room.minY + room.sizeY, ttWall);
        }

        for (let y = 0; y < room.sizeY + 2; ++y) {
            level.set(room.minX - 1, y + room.minY - 1, ttWall);
            level.set(room.minX + room.sizeX, y + room.minY - 1, ttWall);
        }
    }

    // Plot corridors into grid

    for (let roomY = 0; roomY < numCellsY; ++roomY) {
        for (let roomX = 0; roomX < (numCellsX - 1); ++roomX) {
            const roomIndex0 = roomY * numCellsX + roomX;
            const roomIndex1 = roomIndex0 + 1;

            const edge = edges.find(edge => edge.edge[0] === roomIndex0 && edge.edge[1] === roomIndex1);
            if (edge === undefined) {
                continue;
            }

            const edgeCorridorWidth = corridorWidth;

            const room0 = rooms[roomIndex0];
            const room1 = rooms[roomIndex1];

            const xMin = room0.minX + room0.sizeX;
            const xMax = room1.minX;
//            const xMid = randomInRange(xMax - (xMin + 1 + edgeCorridorWidth)) + xMin + 1;
            const xMid = Math.floor((xMax - (xMin + 1 + edgeCorridorWidth)) / 2) + xMin + 1;

            const yMinIntersect = Math.max(room0.minY, room1.minY) + 1;
            const yMaxIntersect = Math.min(room0.minY + room0.sizeY, room1.minY + room1.sizeY) - 1;
            const yRangeIntersect = yMaxIntersect - yMinIntersect;

            let yMinLeft, yMinRight;
            if (yRangeIntersect >= edgeCorridorWidth) {
                yMinLeft = yMinRight = randomInRange(1 + yRangeIntersect - edgeCorridorWidth) + yMinIntersect;
            } else {
                yMinLeft = Math.floor((room0.sizeY - edgeCorridorWidth) / 2) + room0.minY;
                yMinRight = Math.floor((room1.sizeY - edgeCorridorWidth) / 2) + room1.minY;
                // yMinLeft = randomInRange(room0.sizeY - (1 + edgeCorridorWidth)) + room0.minY + 1;
                // yMinRight = randomInRange(room1.sizeY - (1 + edgeCorridorWidth)) + room1.minY + 1;
            }

            for (let x = xMin; x < xMid; ++x) {
                for (let y = 0; y < edgeCorridorWidth; ++y) {
                    level.set(x, yMinLeft + y, ttHall);
                }
            }

            for (let x = xMid + edgeCorridorWidth; x < xMax; ++x) {
                for (let y = 0; y < edgeCorridorWidth; ++y) {
                    level.set(x, yMinRight + y, ttHall);
                }
            }

            const yMin = Math.min(yMinLeft, yMinRight);
            const yMax = Math.max(yMinLeft, yMinRight);
            for (let y = yMin; y < yMax + edgeCorridorWidth; ++y) {
                for (let x = 0; x < edgeCorridorWidth; ++x) {
                    level.set(xMid + x, y, ttHall);
                }
            }
        }
    }

    for (let roomY = 0; roomY < (numCellsY - 1); ++roomY) {
        for (let roomX = 0; roomX < numCellsX; ++roomX) {
            const roomIndex0 = roomY * numCellsX + roomX;
            const roomIndex1 = roomIndex0 + numCellsX;

            const edge = edges.find(edge => edge.edge[0] === roomIndex0 && edge.edge[1] === roomIndex1);
            if (edge === undefined) {
                continue;
            }

            const edgeCorridorWidth = corridorWidth;

            const room0 = rooms[roomIndex0];
            const room1 = rooms[roomIndex1];

            const xMinIntersect = Math.max(room0.minX, room1.minX) + 1;
            const xMaxIntersect = Math.min(room0.minX + room0.sizeX, room1.minX + room1.sizeX) - 1;
            const xRangeIntersect = xMaxIntersect - xMinIntersect;

            let xMinLower, xMinUpper;
            if (xRangeIntersect >= edgeCorridorWidth) {
                xMinLower = xMinUpper = randomInRange(1 + xRangeIntersect - edgeCorridorWidth) + xMinIntersect;
            } else {
                xMinLower = Math.floor((room0.sizeX - edgeCorridorWidth) / 2) + room0.minX;
                xMinUpper = Math.floor((room1.sizeX - edgeCorridorWidth) / 2) + room1.minX;
                // xMinLower = randomInRange(room0.sizeX - (1 + edgeCorridorWidth)) + room0.minX + 1;
                // xMinUpper = randomInRange(room1.sizeX - (1 + edgeCorridorWidth)) + room1.minX + 1;
            }

            const yMin = room0.minY + room0.sizeY;
            const yMax = room1.minY;
//            const yMid = randomInRange(yMax - (yMin + 1 + edgeCorridorWidth)) + yMin + 1;
            const yMid = Math.floor((yMax - (yMin + 1 + edgeCorridorWidth)) / 2) + yMin + 1;

            for (let y = yMin; y < yMid; ++y) {
                for (let x = 0; x < edgeCorridorWidth; ++x) {
                    level.set(xMinLower + x, y, ttHall);
                }
            }

            for (let y = yMid + edgeCorridorWidth; y < yMax; ++y) {
                for (let x = 0; x < edgeCorridorWidth; ++x) {
                    level.set(xMinUpper + x, y, ttHall);
                }
            }

            const xMin = Math.min(xMinLower, xMinUpper);
            const xMax = Math.max(xMinLower, xMinUpper);
            for (let x = xMin; x < xMax + edgeCorridorWidth; ++x) {
                for (let y = 0; y < edgeCorridorWidth; ++y) {
                    level.set(x, yMid + y, ttHall);
                }
            }
        }
    }

    // Convert to colored squares.

    const roomColor = 0xff808080;
    const hallColor = 0xff404040;
    const wallColor = 0xff0055aa;

    const squares = [];
    for (let y = 0; y < level.sizeY; ++y) {
        for (let x = 0; x < level.sizeX; ++x) {
            const type = level.get(x, y);
            if (type == ttRoom) {
                squares.push({x: x, y: y, color: roomColor});
            } else if (type == ttHall) {
                squares.push({x: x, y: y, color: hallColor});
            } else if (type == ttWall) {
                squares.push({x: x, y: y, color: wallColor});
            }
        }
    }

    // Convert squares to triangles

    const numVertices = squares.length * 6;
    const bytesPerVertex = 12;

    const vertexData = new ArrayBuffer(numVertices * bytesPerVertex);
    const vertexDataAsFloat32 = new Float32Array(vertexData);
    const vertexDataAsUint32 = new Uint32Array(vertexData);

    const scale = 1 / 80;

    for (let i = 0; i < squares.length; ++i) {
        const j = 18 * i;
        const color = squares[i].color;
        const x0 = squares[i].x * scale;
        const y0 = squares[i].y * scale;
        const x1 = x0 + scale;
        const y1 = y0 + scale;

        vertexDataAsFloat32[j+0] = x0;
        vertexDataAsFloat32[j+1] = y0;
        vertexDataAsUint32[j+2] = color;

        vertexDataAsFloat32[j+3] = x1;
        vertexDataAsFloat32[j+4] = y0;
        vertexDataAsUint32[j+5] = color;

        vertexDataAsFloat32[j+6] = x0;
        vertexDataAsFloat32[j+7] = y1;
        vertexDataAsUint32[j+8] = color;

        vertexDataAsFloat32[j+9] = x0;
        vertexDataAsFloat32[j+10] = y1;
        vertexDataAsUint32[j+11] = color;

        vertexDataAsFloat32[j+12] = x1;
        vertexDataAsFloat32[j+13] = y0;
        vertexDataAsUint32[j+14] = color;

        vertexDataAsFloat32[j+15] = x1;
        vertexDataAsFloat32[j+16] = y1;
        vertexDataAsUint32[j+17] = color;
    }

    return {
        worldSizeX: level.sizeX * scale,
        worldSizeY: level.sizeY * scale,
        vertexData: vertexData,
    };
}

function shuffleArray(array) {
    for (let i = array.length - 1; i > 0; --i) {
        let j = randomInRange(i + 1);
        let temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}