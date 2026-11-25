#!/bin/bash
# ---------------------------------------------------------------------------
# ProPhosLinker container entrypoint
# 1. Start Neo4j in the background
# 2. Wait for Neo4j to be ready
# 3. Run the ProPhosLinker pipeline with all forwarded arguments
# 4. Stop Neo4j on exit
# ---------------------------------------------------------------------------
set -e

NEO4J_HOME="${NEO4J_HOME:-/opt/neo4j}"
NEO4J_BIN="${NEO4J_HOME}/bin/neo4j"

cleanup() {
    echo ""
    echo "🛑 Shutting down Neo4j..."
    "${NEO4J_BIN}" stop 2>/dev/null || true
    echo "✅ Neo4j stopped."
}
trap cleanup EXIT INT TERM

echo "🚀 Starting Neo4j..."
"${NEO4J_BIN}" start

# Wait for Neo4j HTTP endpoint to become available (max 60s)
echo "⏳ Waiting for Neo4j to be ready..."
for i in $(seq 1 30); do
    if curl -s -o /dev/null -w '%{http_code}' http://localhost:7474 | grep -q '200'; then
        echo "✅ Neo4j is ready (HTTP 200 after ${i}s)"
        break
    fi
    if [ "$i" -eq 30 ]; then
        echo "⚠️  Neo4j didn't start in 30s — continuing anyway"
    fi
    sleep 2
done

echo ""
echo "📊 Running ProPhosLinker..."
echo ""

# Run the pipeline — pass ALL arguments through
exec python3 /app/cli.py "$@"