# For sequences > 100kb, use chunking
if len(seq) > 100000:
    st.warning("Large sequence detected. Using optimized chunked processing...")
    chunk_size = 50000
    overlap = 1000  # Overlap between chunks to avoid missing boundary motifs
    
    all_results = []
    for chunk_start in range(0, len(seq), chunk_size - overlap):
        chunk_end = min(chunk_start + chunk_size, len(seq))
        chunk_seq = seq[chunk_start:chunk_end]
        
        # Process chunk with offset
        chunk_results = process_sequence(chunk_seq, window_size, min_region_len, percentile_thresh)
        
        # Adjust coordinates
        for result in chunk_results:
            result['Motif_Start'] += chunk_start
            result['Motif_End'] += chunk_start
            result['LowP_Start'] += chunk_start
            result['LowP_End'] += chunk_start
        
        all_results.extend(chunk_results)
    
    return all_results
