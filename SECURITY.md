# Security Policy

## Supported versions

Magesty.jl is a numerical Julia package. Only the latest released minor
version receives fixes; older versions are considered unsupported.

| Version | Supported          |
| ------- | ------------------ |
| 0.1.x   | :white_check_mark: |
| < 0.1   | :x:                |

## Reporting a vulnerability

If you believe you have found a security issue (e.g., a vulnerability in a
dependency that affects users of this package, code execution from untrusted
input files, or similar), please **do not** open a public GitHub issue.

Instead, email the maintainer:

- **T. Tanaka** — `tanaka.t.bj3ac@gmail.com`

Please include:

- A description of the issue and its impact.
- Steps to reproduce, ideally with a minimal example.
- Affected version(s).

You should expect an acknowledgement within a few business days. We will
work with you on a coordinated disclosure timeline before any public
discussion.

## Scope

Magesty.jl processes user-supplied TOML inputs, EMBSET data files, and
XML basis/model files. Issues that allow code execution, denial of service,
or data corruption from these inputs are in scope.

Out of scope: general Julia language vulnerabilities, issues in upstream
dependencies (please report those to the respective projects), and
correctness/numerical bugs (please file as a regular GitHub issue).
