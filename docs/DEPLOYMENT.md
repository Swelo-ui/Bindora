# Bindora Dock Production Deployment Guide

## Overview

This guide covers production deployment of Bindora Dock v2.0 with security hardening, database persistence, and containerization.

---

## Quick Start

### Docker Deployment (Recommended)

```bash
# 1. Clone and configure
git clone https://github.com/Swelo-ui/Bindora.git
cd Bindora
cp .env.example .env

# 2. Start services
docker-compose up -d

# 3. Verify deployment
curl http://localhost:5000/api/health
```

### Manual Deployment

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Initialize database
python -c "from backend.db.database import init_db; from backend.config import DATABASE_URL; init_db(DATABASE_URL)"

# 3. Start production server
gunicorn --config gunicorn_config.py backend.app:app
```

---

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     Docker Container                         │
│  ┌────────────────────────────────────────────────────────┐ │
│  │          Gunicorn WSGI (2 workers, 120s timeout)       │ │
│  └────────────────────────────────────────────────────────┘ │
│                             │                                │
│  ┌────────────────────────────────────────────────────────┐ │
│  │                  Flask Application                      │ │
│  │  • Input Validation Layer                              │ │
│  │  • API Routes (/api/*)                                 │ │
│  │  • Database Session Management                         │ │
│  │  • AutoDock Vina Integration                          │ │
│  └────────────────────────────────────────────────────────┘ │
│                             │                                │
│  ┌────────────────────────────────────────────────────────┐ │
│  │              Persistent Volumes                         │ │
│  │  • bindora-data (database, results)                    │ │
│  │  • bindora-cache (PDB structures, API cache)           │ │
│  └────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
```

---

## Configuration

### Environment Variables

Create `.env` file:

```bash
# Security
BINDORA_DEBUG=False
BINDORA_CORS_ORIGINS=http://localhost:5000,http://127.0.0.1:5000
BINDORA_MAX_CONTENT_LENGTH=33554432

# Server
BINDORA_HOST=0.0.0.0
BINDORA_PORT=5000

# Database
DATABASE_URL=sqlite:///data/bindora.db

# Gunicorn (production)
GUNICORN_WORKERS=2
GUNICORN_TIMEOUT=120
```

### Security Checklist

- [ ] `DEBUG=False` in production
- [ ] CORS origins whitelisted
- [ ] Payload size limits configured (32 MB)
- [ ] Database backup strategy in place
- [ ] Container running as non-root user
- [ ] Health checks configured
- [ ] Log rotation enabled

---

## Database Management

### Schema

**Tables:**
- `docking_sessions` - Completed docking simulations
- `batch_screening_jobs` - Batch virtual screening jobs
- `benchmark_records` - Validation benchmark results

### Backup & Restore

```bash
# Backup
docker exec bindora-dock sqlite3 /app/data/bindora.db .dump > backup.sql

# Restore
docker exec -i bindora-dock sqlite3 /app/data/bindora.db < backup.sql

# Copy database out
docker cp bindora-dock:/app/data/bindora.db ./bindora_backup.db
```

### Session Management API

```bash
# List sessions (paginated)
curl http://localhost:5000/api/sessions/list?limit=50&offset=0

# Get session details
curl http://localhost:5000/api/sessions/<session_id>

# Delete session
curl -X DELETE http://localhost:5000/api/sessions/<session_id>
```

---

## Monitoring

### Health Checks

**Endpoint:**
```bash
curl http://localhost:5000/api/health
```

**Response:**
```json
{
  "status": "healthy",
  "service": "Bindora 3D Drug-Receptor & PK/PD Analyzer",
  "vina_available": true,
  "vina_path": "/usr/local/bin/vina"
}
```

### Container Health

```bash
# Check container status
docker ps

# Inspect health
docker inspect bindora-dock | grep -A 10 Health

# View logs
docker-compose logs -f bindora-app
```

### Metrics

Monitor these endpoints:
- `/api/health` - Service availability
- `/api/sessions/list` - Session count growth
- Container CPU/memory usage via `docker stats`

---

## Troubleshooting

### Common Issues

**1. Container fails to start**

```bash
# Check logs
docker-compose logs bindora-app

# Common causes:
# - Port 5000 already in use
# - Vina binary download failed
# - Database initialization error
```

**2. Database locked errors**

```bash
# Stop all processes
docker-compose down

# Remove lock file
docker exec bindora-dock rm /app/data/bindora.db-journal

# Restart
docker-compose up -d
```

**3. Vina not found**

```bash
# Verify Vina binary
docker exec bindora-dock /usr/local/bin/vina --version

# Rebuild container if missing
docker-compose build --no-cache
```

**4. CORS errors in browser**

Update `.env`:
```bash
CORS_ORIGINS=http://localhost:5000,http://127.0.0.1:5000,https://yourdomain.com
```

Restart:
```bash
docker-compose restart
```

---

## Performance Tuning

### Gunicorn Workers

```bash
# Calculate optimal workers
# workers = (2 × CPU_cores) + 1

# For 4-core machine:
GUNICORN_WORKERS=9
```

### Database Optimization

```sql
-- Create indexes (run once)
CREATE INDEX IF NOT EXISTS idx_docking_sessions_timestamp 
  ON docking_sessions(timestamp DESC);
CREATE INDEX IF NOT EXISTS idx_docking_sessions_pdb_id 
  ON docking_sessions(pdb_id);
```

### Cache Management

```bash
# Clear cache (safe to delete)
docker exec bindora-dock rm -rf /app/data/cache/*

# Cache grows over time, consider:
# - Periodic cleanup cron job
# - Volume size monitoring
```

---

## Scaling

### Horizontal Scaling

Use a load balancer (Nginx/Traefik) with multiple containers:

```yaml
# docker-compose.scale.yml
version: '3.8'
services:
  bindora-app:
    # ... configuration
    deploy:
      replicas: 3
  
  nginx:
    image: nginx:alpine
    ports:
      - "80:80"
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf
    depends_on:
      - bindora-app
```

### Database Migration

For production scale, consider migrating from SQLite to PostgreSQL:

```bash
# Update .env
DATABASE_URL=postgresql://user:password@db-host:5432/bindora

# Requires: pip install psycopg2-binary
```

---

## Security Hardening

### Additional Measures

1. **Reverse Proxy** (Nginx/Traefik)
   - SSL/TLS termination
   - Rate limiting
   - DDoS protection

2. **Secrets Management**
   - Use Docker secrets or vault
   - Rotate API keys regularly
   - Never commit .env files

3. **Network Isolation**
   - Run containers in isolated network
   - Firewall rules for port 5000
   - VPN access for admin endpoints

4. **Audit Logging**
   - Log all API requests
   - Monitor failed authentication attempts
   - Alert on anomalous patterns

---

## Maintenance

### Regular Tasks

**Daily:**
- Monitor container health
- Check disk space
- Review error logs

**Weekly:**
- Database backup
- Clear old cache files
- Update dependencies (security patches)

**Monthly:**
- Review session history growth
- Analyze performance metrics
- Test backup restoration

### Update Procedure

```bash
# 1. Backup database
docker cp bindora-dock:/app/data/bindora.db ./backup_$(date +%Y%m%d).db

# 2. Pull latest code
git pull origin main

# 3. Rebuild and restart
docker-compose down
docker-compose build --no-cache
docker-compose up -d

# 4. Verify health
curl http://localhost:5000/api/health
```

---

## Support

For production deployment support:
- Email: sharmaji.pharmatech.info@gmail.com
- GitHub Issues: https://github.com/Swelo-ui/Bindora/issues
- Documentation: README.md, CONTRIBUTING.md

---

## License

Apache 2.0 - See LICENSE file for details.
