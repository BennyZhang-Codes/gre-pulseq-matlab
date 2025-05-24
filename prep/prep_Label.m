function [Label] = prep_Label(Actual, Label, seq)

    % Set PAT scan flag
    % For fast reuse also define an ID for each.

    % Marks GRAPPA line which is used as reference (center of k-space).
    lblSetRefScan            = mr.makeLabel('SET','REF', true ); 
    % Marks GRAPPA line which is used as reference (center of k-space) AND
    % also would be acquired if no reference lines were acquired (e.g., SENSE),
    % i.e., would also be acquired due to the subsampling only.
    lblSetRefAndImaScan      = mr.makeLabel('SET','IMA', true );
    % Negates LabelRef.
    lblResetRefScan          = mr.makeLabel('SET','REF', false);
    % Negates LabelImaAndRef.
    lblResetRefAndImaScan    = mr.makeLabel('SET','IMA', false);

    lblSetRefScan.id         = seq.registerLabelEvent(lblSetRefScan        );
    lblSetRefAndImaScan.id   = seq.registerLabelEvent(lblSetRefAndImaScan  );
    lblResetRefScan.id       = seq.registerLabelEvent(lblResetRefScan      );
    lblResetRefAndImaScan.id = seq.registerLabelEvent(lblResetRefAndImaScan);

    Label.lblSetRefScan            = lblSetRefScan;
    Label.lblSetRefAndImaScan      = lblSetRefAndImaScan;
    Label.lblResetRefScan          = lblResetRefScan;
    Label.lblResetRefAndImaScan    = lblResetRefAndImaScan;
    
    
    % Labels to mark if ADC should be time reversed. Used in bipolar multi-echo acquisitions.
    lblSetRev          = mr.makeLabel('SET', 'REV', true)   ;
    lblSetRev.id       = seq.registerLabelEvent(lblSetRev)  ;
    lblResetRev        = mr.makeLabel('SET', 'REV', false)  ;
    lblResetRev.id     = seq.registerLabelEvent(lblResetRev);

    Label.lblSetRev    = lblSetRev  ;
    Label.lblResetRev  = lblResetRev;
    % Labels to set echo, in case of multi echo acquisition.
    
    % Allocate memory first by creating the last echo label (probably not very important).
    clear LabelEchos ; % clear from previous (possible) run
    % set last echo in array (to initialize structure array)
    LabelEchos(Actual.nTE) = mr.makeLabel('SET', 'ECO', Actual.nTE-1) ; 
    % set all remaining echo labels
    for TECounter = 1:(Actual.nTE-1)
        LabelEchos(TECounter) = mr.makeLabel('SET', 'ECO', TECounter-1) ;
    end
    % Add IDs to all labels (done separately after previous for-loop because
    % once one ID is set the output of mr.makeLabel no longer has the same
    % structure as eLabelEchos(n) (does not contain the 'id' field).
    for TECounter = 1:Actual.nTE
        LabelEchos(TECounter).id = seq.registerLabelEvent(LabelEchos(TECounter)) ;
    end
    Label.lblSetEchos = LabelEchos;
end
