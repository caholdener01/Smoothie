# Copyright (c) 2026 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

"""
Smoothie suppresses certain warnings by default for a better user experience:
- ImplicitModificationWarning: Smoothie intentionally uses views instead of copies
  for memory efficiency when processing large spatial transcriptomics datasets.
- These warnings don't indicate errors - they're expected behavior for the memory-optimized approach.

Users can re-enable warnings with smoothie.enable_warnings() if needed for debugging.
"""

import warnings
import os


def suppress_warnings():
    """
    Suppress specific warnings that Smoothie commonly triggers due to memory-efficient operations.
    
    Only suppresses:
    - AnnData ImplicitModificationWarning (from using views instead of copies)
    - Specific scanpy warnings about working with AnnData views
    - Specific FutureWarning regarding squidpy.pl.spatial_scatter
    
    Does NOT suppress:
    - General numpy/scipy warnings
    - User code warnings
    - Actual errors or problems
    """
    
    # 1. Suppress AnnData ImplicitModificationWarning
    # This occurs when Smoothie intentionally uses views for memory efficiency
    try:
        from anndata import ImplicitModificationWarning
        warnings.filterwarnings('ignore', category=ImplicitModificationWarning)
    except ImportError:
        pass  # If anndata doesn't have this warning class, skip
    
    # 2. Suppress specific scanpy warning about views (only this specific message)
    # This occurs when scanpy receives views that Smoothie intentionally creates
    warnings.filterwarnings(
        'ignore', 
        message='.*Received a view of an AnnData.*',
        category=UserWarning,
        module='scanpy'
    )
    
    # 3. Suppress specific pandas SettingWithCopyWarning when working with AnnData obs/var
    # This can occur during module assignment operations
    try:
        import pandas as pd
        warnings.filterwarnings(
            'ignore',
            category=pd.errors.SettingWithCopyWarning
        )
    except (ImportError, AttributeError):
        pass

    # 4. Suppress specific FutureWarning regarding squidpy replacement
    # Using regex (.*) to handle potential variations in the warning string
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message=".*squidpy.pl.spatial_scatter.*"
    )


def enable_warnings():
    """
    Re-enable all warnings that were suppressed.
    
    This resets ALL warning filters to default, not just Smoothie's.
    """
    warnings.filterwarnings('default')


def quiet_mode():
    """
    Context manager for temporarily suppressing Smoothie warnings.
    
    Use this when you want to suppress warnings for a specific code block only.
    Warnings are automatically restored after the block exits.
    
    Example:
    --------
    >>> with smoothie.quiet_mode():
    ...     # Code here runs without Smoothie warnings
    ...     sm_adata = smoothie.run_gridbased_smoothing(...)
    >>> # Warnings are restored after the block
    """
    class QuietContext:
        def __enter__(self):
            suppress_warnings()
            return self
            
        def __exit__(self, exc_type, exc_val, exc_tb):
            enable_warnings()
            return False
            
    return QuietContext()


# Automatically suppress Smoothie-specific warnings on import by default
# Memory-efficient operations intentionally use views instead of copies
# Users can re-enable with smoothie.enable_warnings() if needed
if os.environ.get('SMOOTHIE_SUPPRESS_WARNINGS', 'true').lower() == 'true':
    suppress_warnings()
