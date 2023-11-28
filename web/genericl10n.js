/* Copyright 2017 Mozilla Foundation
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

class GenericL10n {
  constructor(lang) {
    this._lang = lang;
    this._ready = Promise.resolve();
  }

  async getLanguage() {
    return this._lang;
  }

  async getDirection() {
    // http://www.w3.org/International/questions/qa-scripts
    // Arabic, Hebrew, Farsi, Pashto, Urdu
    const rtlList = ["ar", "he", "fa", "ps", "ur"];
    const shortCode = this._lang.split("-", 1)[0];
    return rtlList.includes(shortCode) ? "rtl" : "ltr";
  }

  async get(key, args, str) {
    if (args && str && str.includes("{")) {
      const reArgs = /\{\{\s*(.+?)\s*\}\}/g;
      str = str.replace(reArgs, function (matched_text, arg) {
        if (arg in args) {
          return args[arg];
        }
        console.warn("argument {{" + arg + "}} for #" + key + " is undefined.");
        return matched_text;
      });
    }

    return str || "";
  }

  async translate() {}
}

export { GenericL10n };
