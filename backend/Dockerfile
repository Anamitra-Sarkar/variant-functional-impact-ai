# syntax=docker/dockerfile:1
# Multi-stage Dockerfile for Variant Impact backend (FastAPI)
# Respects fail-closed release gate: no env vars baked in; defaults to gate CLOSED.

# ---- Builder stage: install Python deps ----
FROM python:3.11-slim AS builder
WORKDIR /app

# System deps for biopython / building (minimal)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# ---- Runtime stage: slim image with installed deps ----
FROM python:3.11-slim AS runtime
WORKDIR /app

# Create non-root user for security
RUN useradd -m appuser

# Copy installed Python packages from builder
COPY --from=builder /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

# Copy application code
# Build context is repo root (see docker-compose.yml), so paths are relative to root
COPY backend/ ./backend/
COPY data_pipeline/ ./data_pipeline/
COPY requirements.txt ./requirements.txt

# Ensure artifacts dir exists (honest empty state when gate closed)
RUN mkdir -p artifacts && chown -R appuser:appuser /app

USER appuser

EXPOSE 8000

# Release gate is fail-closed by default: do NOT set MODEL_RELEASE_APPROVED here.
# It must be set explicitly via environment (docker-compose.yml or deployment).
# ENV MODEL_RELEASE_APPROVED=false  # implicit default is unset -> closed

# Health check uses the always-available /health endpoint
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import urllib.request, sys; sys.exit(0 if urllib.request.urlopen('http://localhost:8000/health').getcode()==200 else 1)" || exit 1

CMD ["uvicorn", "backend.app:app", "--host", "0.0.0.0", "--port", "8000"]
