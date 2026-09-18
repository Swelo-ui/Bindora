"""
Gunicorn WSGI server configuration for production deployment
"""

import os
import multiprocessing

# Server socket
bind = f"{os.getenv('BINDORA_HOST', '0.0.0.0')}:{os.getenv('BINDORA_PORT', '5000')}"

# Worker processes
workers = int(os.getenv('GUNICORN_WORKERS', 2))
worker_class = "sync"
timeout = int(os.getenv('GUNICORN_TIMEOUT', 120))

# Restart workers after processing this many requests (prevent memory leaks)
max_requests = 1000
max_requests_jitter = 50

# Logging
accesslog = "-"  # stdout
errorlog = "-"   # stderr
loglevel = os.getenv('GUNICORN_LOG_LEVEL', 'info')

# Process naming
proc_name = "bindora"

# Server mechanics
daemon = False
pidfile = None
umask = 0
user = None
group = None
tmp_upload_dir = None

# Security
limit_request_line = 4096
limit_request_fields = 100
limit_request_field_size = 8190
