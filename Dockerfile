# Multi-stage Dockerfile for Bindora Dock production deployment
# Python 3.11 slim base with AutoDock Vina 1.2.5 Linux binary

FROM python:3.11-slim AS base

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Download and install AutoDock Vina 1.2.5 Linux binary
FROM base AS vina-downloader
WORKDIR /tmp
RUN wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 \
    && mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina \
    && chmod +x /usr/local/bin/vina \
    && /usr/local/bin/vina --version

# Python dependencies stage
FROM base AS dependencies
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Final application stage
FROM dependencies AS final

# Copy Vina binary from downloader stage
COPY --from=vina-downloader /usr/local/bin/vina /usr/local/bin/vina

# Create non-root user
RUN useradd -m -u 1000 -s /bin/bash bindora \
    && mkdir -p /app/data /app/data/cache /app/data/benchmarks \
    && chown -R bindora:bindora /app

# Copy application code
COPY --chown=bindora:bindora . /app

# Switch to non-root user
USER bindora

# Set working directory
WORKDIR /app

# Expose port
EXPOSE 5000

# Health check
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD wget --no-verbose --tries=1 --spider http://localhost:5000/api/health || exit 1

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    VINA_EXE=/usr/local/bin/vina \
    BINDORA_HOST=0.0.0.0 \
    BINDORA_PORT=5000

# Start Gunicorn WSGI server
CMD ["gunicorn", "--config", "gunicorn_config.py", "backend.app:app"]
