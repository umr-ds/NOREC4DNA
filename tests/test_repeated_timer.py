#!/usr/bin/python
# -*- coding: latin-1 -*-
"""
Coverage tests for RepeatedTimer module.
"""
import time
import pytest

from norec4dna.helper.RepeatedTimer import RepeatedTimer


class TestRepeatedTimer:
    """Test RepeatedTimer class - aims for 100% coverage"""
    
    def test_repeated_timer_basic(self):
        """Test basic RepeatedTimer functionality"""
        call_count = [0]
        
        def increment():
            call_count[0] += 1
        
        # Create timer with 0.1 second interval
        rt = RepeatedTimer(0.1, increment)
        
        # Let it run for a bit
        time.sleep(0.25)
        
        # Stop the timer
        rt.stop()
        
        # Should have been called at least twice
        assert call_count[0] >= 2
    
    def test_repeated_timer_with_args(self):
        """Test RepeatedTimer with arguments"""
        results = []
        
        def add_value(x, y):
            results.append(x + y)
        
        rt = RepeatedTimer(0.1, add_value, 5, 3)
        time.sleep(0.15)
        rt.stop()
        
        assert len(results) >= 1
        assert all(r == 8 for r in results)
    
    def test_repeated_timer_with_kwargs(self):
        """Test RepeatedTimer with keyword arguments"""
        results = []
        
        def add_named(a=0, b=0):
            results.append(a + b)
        
        rt = RepeatedTimer(0.1, add_named, a=10, b=20)
        time.sleep(0.15)
        rt.stop()
        
        assert len(results) >= 1
        assert all(r == 30 for r in results)
    
    def test_repeated_timer_start_stop(self):
        """Test RepeatedTimer start and stop methods"""
        call_count = [0]
        
        def increment():
            call_count[0] += 1
        
        rt = RepeatedTimer(0.1, increment)
        assert rt.is_running
        
        rt.stop()
        assert not rt.is_running
        
        # Restart
        rt.start()
        assert rt.is_running
        
        time.sleep(0.15)
        rt.stop()
        assert call_count[0] >= 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
