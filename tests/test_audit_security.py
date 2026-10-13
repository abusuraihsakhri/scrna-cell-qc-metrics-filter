"""Security regression tests for the HMAC-backed audit record."""

from agents.base import AuditTrail


def test_signature_integrity_detects_tampering():
    log = AuditTrail(secret_key="local-test-key")
    entry = log.log("worker", "qc", "TEST", {"status": "ok"})
    assert log.verify_integrity()
    assert entry["current_hash"]
    log.logs[0]["event_type"] = "FORGED"
    assert not log.verify_integrity()


def test_audit_get_trail_is_defensive_copy():
    log = AuditTrail(secret_key="local-test-key")
    log.log("worker", "qc", "TEST", {"status": "ok"})
    retrieved = log.get_trail()
    retrieved[0]["actor"] = "forged"
    retrieved.clear()
    assert log.verify_integrity()
    assert len(log.get_trail()) == 1
