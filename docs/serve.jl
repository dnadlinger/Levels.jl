# Serve documentation locally, rebuilding it on changes.
#
# For the change detection, skip the generated Literate.jl output pages to avoid loops.

using LiveServer
servedocs(; include_dirs=["src", "examples"], skip_dirs=["docs/src/30-examples"])
