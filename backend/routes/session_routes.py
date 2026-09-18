"""
Session history management API routes
"""

from flask import Blueprint, jsonify, request
from backend.db.database import get_db
from backend.db.models import DockingSession
from sqlalchemy import desc

session_bp = Blueprint('sessions', __name__, url_prefix='/api/sessions')


@session_bp.route('/list', methods=['GET'])
def list_sessions():
    """
    Get paginated list of docking sessions
    
    Query params:
    - limit: max records to return (default: 50)
    - offset: pagination offset (default: 0)
    - pdb_id: filter by PDB ID (optional)
    """
    try:
        limit = int(request.args.get('limit', 50))
        offset = int(request.args.get('offset', 0))
        pdb_id_filter = request.args.get('pdb_id', '').strip()
        
        # Validate limits
        limit = min(max(1, limit), 200)  # Between 1 and 200
        offset = max(0, offset)
        
        with get_db() as db:
            query = db.query(DockingSession).order_by(desc(DockingSession.timestamp))
            
            # Apply filters
            if pdb_id_filter:
                query = query.filter(DockingSession.pdb_id == pdb_id_filter.upper())
            
            # Get total count
            total = query.count()
            
            # Apply pagination
            sessions = query.limit(limit).offset(offset).all()
            
            return jsonify({
                'sessions': [s.to_dict() for s in sessions],
                'total': total,
                'limit': limit,
                'offset': offset
            })
    except Exception as e:
        return jsonify({'error': f'Failed to list sessions: {str(e)}'}), 500


@session_bp.route('/<session_id>', methods=['GET'])
def get_session(session_id):
    """Get detailed session information by ID"""
    try:
        with get_db() as db:
            session = db.query(DockingSession).filter_by(id=session_id).first()
            
            if not session:
                return jsonify({'error': 'Session not found'}), 404
            
            return jsonify(session.to_dict())
    except Exception as e:
        return jsonify({'error': f'Failed to retrieve session: {str(e)}'}), 500


@session_bp.route('/<session_id>', methods=['DELETE'])
def delete_session(session_id):
    """Delete a session by ID"""
    try:
        with get_db() as db:
            session = db.query(DockingSession).filter_by(id=session_id).first()
            
            if not session:
                return jsonify({'error': 'Session not found'}), 404
            
            db.delete(session)
            db.commit()
            
            return '', 204
    except Exception as e:
        return jsonify({'error': f'Failed to delete session: {str(e)}'}), 500


@session_bp.route('/recent', methods=['GET'])
def get_recent_sessions():
    """Get 10 most recent docking sessions"""
    try:
        with get_db() as db:
            sessions = db.query(DockingSession)\
                .order_by(desc(DockingSession.timestamp))\
                .limit(10)\
                .all()
            
            return jsonify({
                'sessions': [s.to_dict() for s in sessions]
            })
    except Exception as e:
        return jsonify({'error': f'Failed to retrieve recent sessions: {str(e)}'}), 500
