"""
	ProgressReporting

Internal helper that wraps `ProgressMeter.Progress` with a TTY-aware
fallback. On an interactive terminal, an `@threads` loop reports
per-iteration progress as a live bar; on a non-TTY stream (log file,
`tee`, CI) the same call site emits a single `"<desc>... Done (X.XX sec)."`
line instead, so build logs remain grep-friendly without an embedded
carriage-return history.
"""
module ProgressReporting

using Printf
using ProgressMeter: Progress, next!, finish!

export with_progress, tick!

"""
	ProgressHandle

Thin wrapper passed to the body of [`with_progress`](@ref). The
underlying `prog` is either a live `ProgressMeter.Progress` (TTY mode)
or `nothing` (non-TTY or disabled). [`tick!`](@ref) does the right thing
in both cases.
"""
struct ProgressHandle
	prog::Union{Progress, Nothing}
end

"""
	tick!(h::ProgressHandle)

Advance the bar by one step. No-op when the handle has no underlying
bar. Safe to call from inside `@threads`; the underlying
`ProgressMeter.next!` takes an internal lock.
"""
@inline tick!(h::ProgressHandle) = (h.prog === nothing || next!(h.prog); nothing)

"""
	with_progress(body, n, desc; verbosity=true)

Run `body(handle)` with progress reporting around an `n`-iteration
block. The block is expected to call `tick!(handle)` once per
iteration.

`verbosity == false` makes the call transparent: `body` runs with a
no-op handle and nothing is printed. Otherwise output mode is chosen
by `isa(stderr, Base.TTY)`:

- TTY: a live `Progress` bar with description `desc`.
- Non-TTY: one `"<desc>... "` line before the block and
  `"Done (X.XX sec)."` after — no bar, no `\\r` history in the log.

Returns whatever `body` returns.
"""
function with_progress(
	body,
	n::Integer,
	desc::AbstractString;
	verbosity::Bool = true,
)
	if !verbosity
		return body(ProgressHandle(nothing))
	end
	if isa(stderr, Base.TTY)
		prog = Progress(n; desc = desc * " ", barlen = 40)
		try
			return body(ProgressHandle(prog))
		finally
			finish!(prog)
		end
	else
		print(stderr, desc * "... ")
		flush(stderr)
		t0 = time_ns()
		try
			return body(ProgressHandle(nothing))
		finally
			elapsed = (time_ns() - t0) / 1e9
			println(stderr, @sprintf("Done (%.2f sec).", elapsed))
		end
	end
end

end # module
