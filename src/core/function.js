/* Copyright 2012 Mozilla Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import { Dict, isDict, isStream, Ref } from "./primitives.js";
import {
  FormatError,
  info,
  isBool,
  IsEvalSupportedCached,
  unreachable,
} from "../shared/util.js";
import { PostScriptLexer, PostScriptParser } from "./ps_parser.js";
import { LocalFunctionCache } from "./image_utils.js";

class PDFFunctionFactory {
  constructor({ xref, isEvalSupported = true }) {
    this.xref = xref;
    this.isEvalSupported = false;
    this._localFunctionCache = null; // Initialized lazily.
  }

  create(fn) {
    const cachedFunction = this.getCached(fn);
    if (cachedFunction) {
      return cachedFunction;
    }
    const parsedFunction = PDFFunction.parse({
      xref: this.xref,
      isEvalSupported: this.isEvalSupported,
      fn: fn instanceof Ref ? this.xref.fetch(fn) : fn,
    });

    // Attempt to cache the parsed Function, by reference.
    this._cache(fn, parsedFunction);

    return parsedFunction;
  }

  createFromArray(fnObj) {
    const cachedFunction = this.getCached(fnObj);
    if (cachedFunction) {
      return cachedFunction;
    }
    const parsedFunction = PDFFunction.parseArray({
      xref: this.xref,
      isEvalSupported: this.isEvalSupported,
      fnObj: fnObj instanceof Ref ? this.xref.fetch(fnObj) : fnObj,
    });

    // Attempt to cache the parsed Function, by reference.
    this._cache(fnObj, parsedFunction);

    return parsedFunction;
  }

  getCached(cacheKey) {
    let fnRef;
    if (cacheKey instanceof Ref) {
      fnRef = cacheKey;
    } else if (cacheKey instanceof Dict) {
      fnRef = cacheKey.objId;
    } else if (isStream(cacheKey)) {
      fnRef = cacheKey.dict && cacheKey.dict.objId;
    }
    if (fnRef) {
      if (!this._localFunctionCache) {
        this._localFunctionCache = new LocalFunctionCache();
      }
      const localFunction = this._localFunctionCache.getByRef(fnRef);
      if (localFunction) {
        return localFunction;
      }
    }
    return null;
  }

  /**
   * @private
   */
  _cache(cacheKey, parsedFunction) {
    if (!parsedFunction) {
      throw new Error(
        'PDFFunctionFactory._cache - expected "parsedFunction" argument.'
      );
    }
    let fnRef;
    if (cacheKey instanceof Ref) {
      fnRef = cacheKey;
    } else if (cacheKey instanceof Dict) {
      fnRef = cacheKey.objId;
    } else if (isStream(cacheKey)) {
      fnRef = cacheKey.dict && cacheKey.dict.objId;
    }
    if (fnRef) {
      if (!this._localFunctionCache) {
        this._localFunctionCache = new LocalFunctionCache();
      }
      this._localFunctionCache.set(/* name = */ null, fnRef, parsedFunction);
    }
  }
}

function toNumberArray(arr) {
  if (!Array.isArray(arr)) {
    return null;
  }
  const length = arr.length;
  for (let i = 0; i < length; i++) {
    if (typeof arr[i] !== "number") {
      // Non-number is found -- convert all items to numbers.
      const result = new Array(length);
      for (let j = 0; j < length; j++) {
        result[j] = +arr[j];
      }
      return result;
    }
  }
  return arr;
}

var PDFFunction = (function PDFFunctionClosure() {
  const CONSTRUCT_SAMPLED = 0;
  const CONSTRUCT_INTERPOLATED = 2;
  const CONSTRUCT_STICHED = 3;
  const CONSTRUCT_POSTSCRIPT = 4;

  return {
    getSampleArray(size, outputSize, bps, stream) {
      var i, ii;
      var length = 1;
      for (i = 0, ii = size.length; i < ii; i++) {
        length *= size[i];
      }
      length *= outputSize;

      var array = new Array(length);
      var codeSize = 0;
      var codeBuf = 0;
      // 32 is a valid bps so shifting won't work
      var sampleMul = 1.0 / (2.0 ** bps - 1);

      var strBytes = stream.getBytes((length * bps + 7) / 8);
      var strIdx = 0;
      for (i = 0; i < length; i++) {
        while (codeSize < bps) {
          codeBuf <<= 8;
          codeBuf |= strBytes[strIdx++];
          codeSize += 8;
        }
        codeSize -= bps;
        array[i] = (codeBuf >> codeSize) * sampleMul;
        codeBuf &= (1 << codeSize) - 1;
      }
      return array;
    },

    getIR({ xref, isEvalSupported, fn }) {
      var dict = fn.dict;
      if (!dict) {
        dict = fn;
      }

      var types = [
        this.constructSampled,
        null,
        this.constructInterpolated,
        this.constructStiched,
        this.constructPostScript,
      ];

      var typeNum = dict.get("FunctionType");
      var typeFn = types[typeNum];
      if (!typeFn) {
        throw new FormatError("Unknown type of function");
      }

      return typeFn.call(this, { xref, isEvalSupported, fn, dict });
    },

    fromIR({ xref, isEvalSupported, IR }) {
      var type = IR[0];
      switch (type) {
        case CONSTRUCT_SAMPLED:
          return this.constructSampledFromIR({ xref, isEvalSupported, IR });
        case CONSTRUCT_INTERPOLATED:
          return this.constructInterpolatedFromIR({
            xref,
            isEvalSupported,
            IR,
          });
        case CONSTRUCT_STICHED:
          return this.constructStichedFromIR({ xref, isEvalSupported, IR });
        // case CONSTRUCT_POSTSCRIPT:
        default:
          return this.constructPostScriptFromIR({ xref, isEvalSupported, IR });
      }
    },

    parse({ xref, isEvalSupported, fn }) {
      const IR = this.getIR({ xref, isEvalSupported, fn });
      return this.fromIR({ xref, isEvalSupported, IR });
    },

    parseArray({ xref, isEvalSupported, fnObj }) {
      if (!Array.isArray(fnObj)) {
        // not an array -- parsing as regular function
        return this.parse({ xref, isEvalSupported, fn: fnObj });
      }

      var fnArray = [];
      for (var j = 0, jj = fnObj.length; j < jj; j++) {
        fnArray.push(
          this.parse({ xref, isEvalSupported, fn: xref.fetchIfRef(fnObj[j]) })
        );
      }
      return function (src, srcOffset, dest, destOffset) {
        for (var i = 0, ii = fnArray.length; i < ii; i++) {
          fnArray[i](src, srcOffset, dest, destOffset + i);
        }
      };
    },

    constructSampled({ xref, isEvalSupported, fn, dict }) {
      function toMultiArray(arr) {
        var inputLength = arr.length;
        var out = [];
        var index = 0;
        for (var i = 0; i < inputLength; i += 2) {
          out[index] = [arr[i], arr[i + 1]];
          ++index;
        }
        return out;
      }
      var domain = toNumberArray(dict.getArray("Domain"));
      var range = toNumberArray(dict.getArray("Range"));

      if (!domain || !range) {
        throw new FormatError("No domain or range");
      }

      var inputSize = domain.length / 2;
      var outputSize = range.length / 2;

      domain = toMultiArray(domain);
      range = toMultiArray(range);

      var size = toNumberArray(dict.getArray("Size"));
      var bps = dict.get("BitsPerSample");
      var order = dict.get("Order") || 1;
      if (order !== 1) {
        // No description how cubic spline interpolation works in PDF32000:2008
        // As in poppler, ignoring order, linear interpolation may work as good
        info("No support for cubic spline interpolation: " + order);
      }

      var encode = toNumberArray(dict.getArray("Encode"));
      if (!encode) {
        encode = [];
        for (var i = 0; i < inputSize; ++i) {
          encode.push([0, size[i] - 1]);
        }
      } else {
        encode = toMultiArray(encode);
      }

      var decode = toNumberArray(dict.getArray("Decode"));
      if (!decode) {
        decode = range;
      } else {
        decode = toMultiArray(decode);
      }

      var samples = this.getSampleArray(size, outputSize, bps, fn);

      return [
        CONSTRUCT_SAMPLED,
        inputSize,
        domain,
        encode,
        decode,
        samples,
        size,
        outputSize,
        2 ** bps - 1,
        range,
      ];
    },

    constructSampledFromIR({ xref, isEvalSupported, IR }) {
      // See chapter 3, page 109 of the PDF reference
      function interpolate(x, xmin, xmax, ymin, ymax) {
        return ymin + (x - xmin) * ((ymax - ymin) / (xmax - xmin));
      }

      return function constructSampledFromIRResult(
        src,
        srcOffset,
        dest,
        destOffset
      ) {
        // See chapter 3, page 110 of the PDF reference.
        var m = IR[1];
        var domain = IR[2];
        var encode = IR[3];
        var decode = IR[4];
        var samples = IR[5];
        var size = IR[6];
        var n = IR[7];
        // var mask = IR[8];
        var range = IR[9];

        // Building the cube vertices: its part and sample index
        // http://rjwagner49.com/Mathematics/Interpolation.pdf
        var cubeVertices = 1 << m;
        var cubeN = new Float64Array(cubeVertices);
        var cubeVertex = new Uint32Array(cubeVertices);
        var i, j;
        for (j = 0; j < cubeVertices; j++) {
          cubeN[j] = 1;
        }

        var k = n,
          pos = 1;
        // Map x_i to y_j for 0 <= i < m using the sampled function.
        for (i = 0; i < m; ++i) {
          // x_i' = min(max(x_i, Domain_2i), Domain_2i+1)
          var domain_2i = domain[i][0];
          var domain_2i_1 = domain[i][1];
          var xi = Math.min(
            Math.max(src[srcOffset + i], domain_2i),
            domain_2i_1
          );

          // e_i = Interpolate(x_i', Domain_2i, Domain_2i+1,
          //                   Encode_2i, Encode_2i+1)
          var e = interpolate(
            xi,
            domain_2i,
            domain_2i_1,
            encode[i][0],
            encode[i][1]
          );

          // e_i' = min(max(e_i, 0), Size_i - 1)
          var size_i = size[i];
          e = Math.min(Math.max(e, 0), size_i - 1);

          // Adjusting the cube: N and vertex sample index
          var e0 = e < size_i - 1 ? Math.floor(e) : e - 1; // e1 = e0 + 1;
          var n0 = e0 + 1 - e; // (e1 - e) / (e1 - e0);
          var n1 = e - e0; // (e - e0) / (e1 - e0);
          var offset0 = e0 * k;
          var offset1 = offset0 + k; // e1 * k
          for (j = 0; j < cubeVertices; j++) {
            if (j & pos) {
              cubeN[j] *= n1;
              cubeVertex[j] += offset1;
            } else {
              cubeN[j] *= n0;
              cubeVertex[j] += offset0;
            }
          }

          k *= size_i;
          pos <<= 1;
        }

        for (j = 0; j < n; ++j) {
          // Sum all cube vertices' samples portions
          var rj = 0;
          for (i = 0; i < cubeVertices; i++) {
            rj += samples[cubeVertex[i] + j] * cubeN[i];
          }

          // r_j' = Interpolate(r_j, 0, 2^BitsPerSample - 1,
          //                    Decode_2j, Decode_2j+1)
          rj = interpolate(rj, 0, 1, decode[j][0], decode[j][1]);

          // y_j = min(max(r_j, range_2j), range_2j+1)
          dest[destOffset + j] = Math.min(
            Math.max(rj, range[j][0]),
            range[j][1]
          );
        }
      };
    },

    constructInterpolated({ xref, isEvalSupported, fn, dict }) {
      var c0 = toNumberArray(dict.getArray("C0")) || [0];
      var c1 = toNumberArray(dict.getArray("C1")) || [1];
      var n = dict.get("N");

      var length = c0.length;
      var diff = [];
      for (var i = 0; i < length; ++i) {
        diff.push(c1[i] - c0[i]);
      }

      return [CONSTRUCT_INTERPOLATED, c0, diff, n];
    },

    constructInterpolatedFromIR({ xref, isEvalSupported, IR }) {
      var c0 = IR[1];
      var diff = IR[2];
      var n = IR[3];

      var length = diff.length;

      return function constructInterpolatedFromIRResult(
        src,
        srcOffset,
        dest,
        destOffset
      ) {
        var x = n === 1 ? src[srcOffset] : src[srcOffset] ** n;

        for (var j = 0; j < length; ++j) {
          dest[destOffset + j] = c0[j] + x * diff[j];
        }
      };
    },

    constructStiched({ xref, isEvalSupported, fn, dict }) {
      var domain = toNumberArray(dict.getArray("Domain"));

      if (!domain) {
        throw new FormatError("No domain");
      }

      var inputSize = domain.length / 2;
      if (inputSize !== 1) {
        throw new FormatError("Bad domain for stiched function");
      }

      var fnRefs = dict.get("Functions");
      var fns = [];
      for (var i = 0, ii = fnRefs.length; i < ii; ++i) {
        fns.push(
          this.parse({ xref, isEvalSupported, fn: xref.fetchIfRef(fnRefs[i]) })
        );
      }

      var bounds = toNumberArray(dict.getArray("Bounds"));
      var encode = toNumberArray(dict.getArray("Encode"));

      return [CONSTRUCT_STICHED, domain, bounds, encode, fns];
    },

    constructStichedFromIR({ xref, isEvalSupported, IR }) {
      var domain = IR[1];
      var bounds = IR[2];
      var encode = IR[3];
      var fns = IR[4];
      var tmpBuf = new Float32Array(1);

      return function constructStichedFromIRResult(
        src,
        srcOffset,
        dest,
        destOffset
      ) {
        var clip = function constructStichedFromIRClip(v, min, max) {
          if (v > max) {
            v = max;
          } else if (v < min) {
            v = min;
          }
          return v;
        };

        // clip to domain
        var v = clip(src[srcOffset], domain[0], domain[1]);
        // calculate which bound the value is in
        for (var i = 0, ii = bounds.length; i < ii; ++i) {
          if (v < bounds[i]) {
            break;
          }
        }

        // encode value into domain of function
        var dmin = domain[0];
        if (i > 0) {
          dmin = bounds[i - 1];
        }
        var dmax = domain[1];
        if (i < bounds.length) {
          dmax = bounds[i];
        }

        var rmin = encode[2 * i];
        var rmax = encode[2 * i + 1];

        // Prevent the value from becoming NaN as a result
        // of division by zero (fixes issue6113.pdf).
        tmpBuf[0] =
          dmin === dmax
            ? rmin
            : rmin + ((v - dmin) * (rmax - rmin)) / (dmax - dmin);

        // call the appropriate function
        fns[i](tmpBuf, 0, dest, destOffset);
      };
    },

    constructPostScript({ xref, isEvalSupported, fn, dict }) {
      var domain = toNumberArray(dict.getArray("Domain"));
      var range = toNumberArray(dict.getArray("Range"));

      if (!domain) {
        throw new FormatError("No domain.");
      }

      if (!range) {
        throw new FormatError("No range.");
      }

      var lexer = new PostScriptLexer(fn);
      var parser = new PostScriptParser(lexer);
      var code = parser.parse();

      return [CONSTRUCT_POSTSCRIPT, domain, range, code];
    },

    constructPostScriptFromIR({ xref, isEvalSupported, IR }) {
      var domain = IR[1];
      var range = IR[2];
      var code = IR[3];

      var numOutputs = range.length >> 1;
      var numInputs = domain.length >> 1;
      var evaluator = new PostScriptEvaluator(code);
      // Cache the values for a big speed up, the cache size is limited though
      // since the number of possible values can be huge from a PS function.
      var cache = Object.create(null);
      // The MAX_CACHE_SIZE is set to ~4x the maximum number of distinct values
      // seen in our tests.
      var MAX_CACHE_SIZE = 2048 * 4;
      var cache_available = MAX_CACHE_SIZE;
      var tmpBuf = new Float32Array(numInputs);

      return function constructPostScriptFromIRResult(
        src,
        srcOffset,
        dest,
        destOffset
      ) {
        var i, value;
        var key = "";
        var input = tmpBuf;
        for (i = 0; i < numInputs; i++) {
          value = src[srcOffset + i];
          input[i] = value;
          key += value + "_";
        }

        var cachedValue = cache[key];
        if (cachedValue !== undefined) {
          dest.set(cachedValue, destOffset);
          return;
        }

        var output = new Float32Array(numOutputs);
        var stack = evaluator.execute(input);
        var stackIndex = stack.length - numOutputs;
        for (i = 0; i < numOutputs; i++) {
          value = stack[stackIndex + i];
          var bound = range[i * 2];
          if (value < bound) {
            value = bound;
          } else {
            bound = range[i * 2 + 1];
            if (value > bound) {
              value = bound;
            }
          }
          output[i] = value;
        }
        if (cache_available > 0) {
          cache_available--;
          cache[key] = output;
        }
        dest.set(output, destOffset);
      };
    },
  };
})();

function isPDFFunction(v) {
  var fnDict;
  if (typeof v !== "object") {
    return false;
  } else if (isDict(v)) {
    fnDict = v;
  } else if (isStream(v)) {
    fnDict = v.dict;
  } else {
    return false;
  }
  return fnDict.has("FunctionType");
}

var PostScriptStack = (function PostScriptStackClosure() {
  var MAX_STACK_SIZE = 100;

  // eslint-disable-next-line no-shadow
  function PostScriptStack(initialStack) {
    this.stack = !initialStack
      ? []
      : Array.prototype.slice.call(initialStack, 0);
  }

  PostScriptStack.prototype = {
    push: function PostScriptStack_push(value) {
      if (this.stack.length >= MAX_STACK_SIZE) {
        throw new Error("PostScript function stack overflow.");
      }
      this.stack.push(value);
    },
    pop: function PostScriptStack_pop() {
      if (this.stack.length <= 0) {
        throw new Error("PostScript function stack underflow.");
      }
      return this.stack.pop();
    },
    copy: function PostScriptStack_copy(n) {
      if (this.stack.length + n >= MAX_STACK_SIZE) {
        throw new Error("PostScript function stack overflow.");
      }
      var stack = this.stack;
      for (var i = stack.length - n, j = n - 1; j >= 0; j--, i++) {
        stack.push(stack[i]);
      }
    },
    index: function PostScriptStack_index(n) {
      this.push(this.stack[this.stack.length - n - 1]);
    },
    // rotate the last n stack elements p times
    roll: function PostScriptStack_roll(n, p) {
      var stack = this.stack;
      var l = stack.length - n;
      var r = stack.length - 1,
        c = l + (p - Math.floor(p / n) * n),
        i,
        j,
        t;
      for (i = l, j = r; i < j; i++, j--) {
        t = stack[i];
        stack[i] = stack[j];
        stack[j] = t;
      }
      for (i = l, j = c - 1; i < j; i++, j--) {
        t = stack[i];
        stack[i] = stack[j];
        stack[j] = t;
      }
      for (i = c, j = r; i < j; i++, j--) {
        t = stack[i];
        stack[i] = stack[j];
        stack[j] = t;
      }
    },
  };
  return PostScriptStack;
})();
var PostScriptEvaluator = (function PostScriptEvaluatorClosure() {
  // eslint-disable-next-line no-shadow
  function PostScriptEvaluator(operators) {
    this.operators = operators;
  }
  PostScriptEvaluator.prototype = {
    execute: function PostScriptEvaluator_execute(initialStack) {
      var stack = new PostScriptStack(initialStack);
      var counter = 0;
      var operators = this.operators;
      var length = operators.length;
      var operator, a, b;
      while (counter < length) {
        operator = operators[counter++];
        if (typeof operator === "number") {
          // Operator is really an operand and should be pushed to the stack.
          stack.push(operator);
          continue;
        }
        switch (operator) {
          // non standard ps operators
          case "jz": // jump if false
            b = stack.pop();
            a = stack.pop();
            if (!a) {
              counter = b;
            }
            break;
          case "j": // jump
            a = stack.pop();
            counter = a;
            break;

          // all ps operators in alphabetical order (excluding if/ifelse)
          case "abs":
            a = stack.pop();
            stack.push(Math.abs(a));
            break;
          case "add":
            b = stack.pop();
            a = stack.pop();
            stack.push(a + b);
            break;
          case "and":
            b = stack.pop();
            a = stack.pop();
            if (isBool(a) && isBool(b)) {
              stack.push(a && b);
            } else {
              stack.push(a & b);
            }
            break;
          case "atan":
            a = stack.pop();
            stack.push(Math.atan(a));
            break;
          case "bitshift":
            b = stack.pop();
            a = stack.pop();
            if (a > 0) {
              stack.push(a << b);
            } else {
              stack.push(a >> b);
            }
            break;
          case "ceiling":
            a = stack.pop();
            stack.push(Math.ceil(a));
            break;
          case "copy":
            a = stack.pop();
            stack.copy(a);
            break;
          case "cos":
            a = stack.pop();
            stack.push(Math.cos(a));
            break;
          case "cvi":
            a = stack.pop() | 0;
            stack.push(a);
            break;
          case "cvr":
            // noop
            break;
          case "div":
            b = stack.pop();
            a = stack.pop();
            stack.push(a / b);
            break;
          case "dup":
            stack.copy(1);
            break;
          case "eq":
            b = stack.pop();
            a = stack.pop();
            stack.push(a === b);
            break;
          case "exch":
            stack.roll(2, 1);
            break;
          case "exp":
            b = stack.pop();
            a = stack.pop();
            stack.push(a ** b);
            break;
          case "false":
            stack.push(false);
            break;
          case "floor":
            a = stack.pop();
            stack.push(Math.floor(a));
            break;
          case "ge":
            b = stack.pop();
            a = stack.pop();
            stack.push(a >= b);
            break;
          case "gt":
            b = stack.pop();
            a = stack.pop();
            stack.push(a > b);
            break;
          case "idiv":
            b = stack.pop();
            a = stack.pop();
            stack.push((a / b) | 0);
            break;
          case "index":
            a = stack.pop();
            stack.index(a);
            break;
          case "le":
            b = stack.pop();
            a = stack.pop();
            stack.push(a <= b);
            break;
          case "ln":
            a = stack.pop();
            stack.push(Math.log(a));
            break;
          case "log":
            a = stack.pop();
            stack.push(Math.log(a) / Math.LN10);
            break;
          case "lt":
            b = stack.pop();
            a = stack.pop();
            stack.push(a < b);
            break;
          case "mod":
            b = stack.pop();
            a = stack.pop();
            stack.push(a % b);
            break;
          case "mul":
            b = stack.pop();
            a = stack.pop();
            stack.push(a * b);
            break;
          case "ne":
            b = stack.pop();
            a = stack.pop();
            stack.push(a !== b);
            break;
          case "neg":
            a = stack.pop();
            stack.push(-a);
            break;
          case "not":
            a = stack.pop();
            if (isBool(a)) {
              stack.push(!a);
            } else {
              stack.push(~a);
            }
            break;
          case "or":
            b = stack.pop();
            a = stack.pop();
            if (isBool(a) && isBool(b)) {
              stack.push(a || b);
            } else {
              stack.push(a | b);
            }
            break;
          case "pop":
            stack.pop();
            break;
          case "roll":
            b = stack.pop();
            a = stack.pop();
            stack.roll(a, b);
            break;
          case "round":
            a = stack.pop();
            stack.push(Math.round(a));
            break;
          case "sin":
            a = stack.pop();
            stack.push(Math.sin(a));
            break;
          case "sqrt":
            a = stack.pop();
            stack.push(Math.sqrt(a));
            break;
          case "sub":
            b = stack.pop();
            a = stack.pop();
            stack.push(a - b);
            break;
          case "true":
            stack.push(true);
            break;
          case "truncate":
            a = stack.pop();
            a = a < 0 ? Math.ceil(a) : Math.floor(a);
            stack.push(a);
            break;
          case "xor":
            b = stack.pop();
            a = stack.pop();
            if (isBool(a) && isBool(b)) {
              stack.push(a !== b);
            } else {
              stack.push(a ^ b);
            }
            break;
          default:
            throw new FormatError(`Unknown operator ${operator}`);
        }
      }
      return stack.stack;
    },
  };
  return PostScriptEvaluator;
})();

export { isPDFFunction, PDFFunctionFactory, PostScriptEvaluator };
