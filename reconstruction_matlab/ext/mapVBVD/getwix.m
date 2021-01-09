classdef getwix < handle
  %% GETWIX -- read syngo VA35, VB17 and VE11 twix files
  % Author: Christian.Labadie@charite.de
  %         3-nov-2009 (VB), 29-jul-2015 (class), 16-dec-2015 (VA/VD)
  %
  % h = getwix(filename)                Open twix file
  % h = getwix('','map') or h.map()     Map twix file in h.pub.map array
  % h = getwix([],'read')               Store all Mdh/Cdh in h.pub.DMA
  % getwix('','view')                   Preview k-space
  %
  % The handle class getwix(filename) avoids repetitive disk access by 
  % sequentially typecasting from an 8 MB buffer one DMA corresponding to 
  % one readout (see fetchReadout method). The implementation as a handle
  % facilitates the reading of large twix files by allowing to alternate
  % reading and processing of DMAs, see Example (reading EPI).
  % The measurement header and ADC channels are respectively stored in the
  % Mdh and Cdh properties. The DMA info mask may be conveniently
  % interpreted with the flag methods such as MDH_REFLECT or
  % MDH_NOISEADJSCAN. Loop counters are converted from the C++ array
  % convention (0 .. n-1) to that of Matlab (1 .. n) using the lin, acq,
  % par, eco, phs, rep, slc, set, seg, ida, idb, idc, idd, ide methods.
  % Help functions added at the end of this file allow to preview k-space
  % as a stack of concatenated channels, see getwix('','view'), to load all
  % DMAs in the pub property, see getwix(filename,'read'), or to store the 
  % entire twix in the array pub.map and the corresponding masks in 
  % pub.evalInfoMask, see the map method or getwix(filename,'map').
  % In VD/VE format, a multipart twix file is opened assuming that the 
  % relevant data is stored in the last part, see setTwixPart method for 
  % choosing another part.
  % In VA format, the meas.out file should be provided and optionally along
  % with the meas.asc (MrProt) and EVA95.tmp (EVA XProt). 
  % The interpretation of XProt is slow, and should be explicitly requested
  % using the getXProt method. DMAs containing SYNCDATA are not interpreted
  % and should be skipped (see MDH_SYNCDATA).
  %
  %% Example: reading a single slice EPI
  %
  %          r = getwix('','open');
  %          while( r.fetchReadout() > 0 && ~h.MDH_ACQEND() )
  %            for c = 1:r.Mdh.UsedChannels
  %              adc = r.Cdh(c).adc;
  %              if r.MDH_REFLECT(), adc = adc(end:-1:1); end
  %              if r.MDH_PHASCOR()
  %                if r.seg == 1
  %                  r.pub.pc.forward(c, r.acq, :) = adc;
  %                else
  %                  r.pub.pc.reflect(c, r.acq, :) = adc;
  %                end
  %              else
  %                r.pub.ksp(c, r.lin,:) = adc;
  %              end
  %            end
  %          end
  %          r.close();
  %
  %% Properties:
  %    File           - File information:
  %                     fullname, path, name, basename, size, permission
  %                     machinefmt, oskip (size of header in bytes)
  %    MrProt         - Sequence protocol parameters
  %    XProt          - XProtocol parameters
  %    Mdh            - Mdh (array)
  %    pub            - Free holder (e.g. to store reconstruction data)
  %
  %% Methods:
  %    getwix         - Constructor
  %    open           - Open twix file
  %    close          - Close file descriptor without deleting object
  %    setTwixPart(n) - Select from a multi-part twix file (default: last)
  %    fetchReadout   - Stores coil elements for one readout
  %    getXProt       - Process XProt header (e.g. YAPS)
  %    setBuffer()    - Typecast from buffer (default)
  %    setBuffer(MB)  - Change buffer size in MB (default 8MB)
  %    setBuffer(0)   - Reads directly from disk
  %    setWaitbar     - Activate progress bar (e.g. for large twix file)
  %    ok             - Informs if twix file is correctly opened
  %    version        - Date of getwix version
  %
  %% Sequence Counters:
  %    lin            - line index                              (Clin + 1)
  %    acq            - acquisition index                       (Cacq + 1)
  %    par            - partition index                         (Cpar + 1)
  %    eco            - echo index                              (Ceco + 1)
  %    phs            - phase index                             (Cphs + 1)
  %    rep            - measurement repeat index                (Crep + 1)
  %    slc            - slice index                             (Cslc + 1)
  %    set            - set index                               (Cset + 1)
  %    seg            - segment index, TSE or phase correction  (Cseg + 1)
  %    ida            - IceDimension a index                    (Cida + 1)
  %    idb            - IceDimension b index                    (Cidb + 1)
  %    idc            - IceDimension c index                    (Cidc + 1)
  %    idd            - IceDimension d index                    (Cidd + 1)
  %    ide            - IceDimension e index                    (Cide + 1)
  %
  %% Measurement Data Header (MDH) info mask:
  %    infoMDH                - Interpretation of info mask as string
  %    MDH_ACQEND             -  last scan
  %    MDH_RTFEEDBACK         -  Realtime feedback scan
  %    MDH_HPFEEDBACK         -  High perfomance feedback scan
  %    MDH_ONLINE             -  processing should be done online
  %    MDH_OFFLINE            -  processing should be done offline
  %    MDH_SYNCDATA           -  readout contains synchroneous data
  %    MDH_LASTSCANINCONCAT   -  Flag for last scan in concatination
  %    MDH_RAWDATACORRECTION  -  Correct the rawadata with the rawdata correction factor
  %    MDH_LASTSCANINMEAS     -  Flag for last scan in measurement
  %    MDH_SCANSCALEFACTOR    -  Flag for scan specific additional scale factor
  %    MDH_2NDHADAMARPULSE    -  2nd RF exitation of HADAMAR
  %    MDH_REFPHASESTABSCAN   -  reference phase stabilization scan
  %    MDH_PHASESTABSCAN      -  phase stabilization scan
  %    MDH_D3FFT              -  execute 3D FFT
  %    MDH_SIGNREV            -  sign reversal
  %    MDH_PHASEFFT           -  execute phase fft
  %    MDH_SWAPPED            -  swapped phase/readout direction
  %    MDH_POSTSHAREDLINE     -  shared line
  %    MDH_PHASCOR            -  phase correction data
  %    MDH_PATREFSCAN         -  additonal scan for PAT reference line/partition
  %    MDH_PATREFANDIMASCAN   -  additonal scan for PAT reference line/partition that is also used as image scan
  %    MDH_REFLECT            -  reflect line
  %    MDH_NOISEADJSCAN       -  noise adjust scan
  %    MDH_SHARENOW           -  all lines are acquired from the actual and previous e.g. phases
  %    MDH_LASTMEASUREDLINE   -  indicates that the current line is the last measured line of all succeeding e.g. phases
  %    MDH_FIRSTSCANINSLICE   -  indicates first scan in slice (needed for time stamps)
  %    MDH_LASTSCANINSLICE    -  indicates  last scan in slice (needed for time stamps)
  %    MDH_TREFFECTIVEBEGIN   -  indicates the begin time stamp for TReff (triggered measurement)
  %    MDH_TREFFECTIVEEND     -  indicates the   end time stamp for TReff (triggered measurement)
  %    MDH_MDS_REF_POSITION           - indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
  %    MDH_SLC_AVERAGED               - indicates avveraged slice for slice partial averaging scheme
  %    MDH_TAGFLAG1                   - adjust scan
  %    MDH_CT_NORMALIZE               - Marks scans used to calculate correction maps for TimCT-Prescan normalize
  %    MDH_SCAN_FIRST                 - Marks the first scan of a particular map
  %    MDH_SCAN_LAST                  - Marks the last scan of a particular map
  %    MDH_FIRST_SCAN_IN_BLADE        - Marks the first line of a blade
  %    MDH_LAST_SCAN_IN_BLADE         - Marks the last line of a blade
  %    MDH_LAST_BLADE_IN_TR           - Set for all lines of the last BLADE in each TR interval
  %    MDH_PACE                       - Distinguishes PACE scans from non PACE scans
  %    MDH_RETRO_LASTPHASE            - Marks the last phase in a heartbeat
  %    MDH_RETRO_ENDOFMEAS            - Marks an ADC at the end of the measurement
  %    MDH_RETRO_REPEATTHISHEARTBEAT  - Repeat the current heartbeat when this bit is found
  %    MDH_RETRO_REPEATPREVHEARTBEAT  - Repeat the previous heartbeat when this bit is found
  %    MDH_RETRO_ABORTSCANNOW         - Just abort everything
  %    MDH_RETRO_LASTHEARTBEAT        - This adc is from the last heartbeat (a dummy)
  %    MDH_RETRO_DUMMYSCAN            - This adc is just a dummy scan, throw it away
  %    MDH_RETRO_ARRDETDISABLED       - Disable all arrhythmia detection when this bit is found
  %    MDH_B1_CONTROLLOOP             - Marks the readout as to be used for B1 Control Loop
  %    MDH_SKIP_ONLINE_PHASCOR        - Marks scans not to be online phase corrected, even if online phase correction is switched on
  %    MDH_SKIP_REGRIDDING            - Marks scans not to be regridded, even if regridding is switched on
  %
  % For more information: <a href="matlab:doc getwix">doc getwix</a>

  % // Properties /////////////////////////////////////////////////////////
  properties(GetAccess = 'public', SetAccess = 'private')
    File          = [];
    MrProt        = [];
    XProt         = [];
    Mdh           = [];
    Cdh           = [];
  end

  properties(GetAccess = 'private', SetAccess = 'private')
    fid           = -1;
    BufferedDMA   = false;
    BufferPos     = 0;
    BufferSiz     = 0;
    Buffer        = [];
    BufferMxByte  = 0;
    WithWaitBar   = -1;
    WaitBar       = [];
    
	syngoVA		  = 1;
	syngoVB		  = 2;
	syngoVD		  = 4;
    
    Version       = 'getwix, 6-Jan-2016 Christian.Labadie@charite.de';
  end  
	
  properties(GetAccess = 'public', SetAccess = 'public')
    pub           = []; % public data field provided for convenience
  end

  % // Methods ////////////////////////////////////////////////////////////
  methods

    function setTwixPart(h, SelectedPart)
    %% select a part from a multi-part twix file
        fname = h.File.fullname;
        close(h);
        open(h, fname, SelectedPart);
    end
      
    function h = getwix(fname, option)
    %% getwix(fname) -- constructor of getwix open TWIX file and protocol
    %
    % See also: open, setBuffer, setWaitbar, fetchReadout, close
    
        h.fid = -1;
        if nargin == 0 || ~exist(fname,'file')
            h.File.path = pwd;
            if nargin==1 && isdir(fname), h.File.path = fname; end
            [fname, h.File.path] = uigetfile({'meas*.dat;meas*.out','VB/VE twix or VA meas.out file'},'Select Twix file', h.File.path);
            if ~ischar(fname), return; end
            if nargin == 0
                opt = menu({fname; ''; 'Option:'}, ...
                              'view (pub.ksp)', ...
                              'read (pub.DMA)', ...
                              'map (pub.map)', ...
                              'open', ...
                              'quit');
                if opt == 5
                    clear
                    return
                elseif opt == 3
                    option = 'map';
                elseif opt == 2
                    option = 'read';
                elseif opt == 1
                    option = 'view';
                end
            end
            fname = fullfile(h.File.path,fname);
        end
        
        open(h, fname, -1);
        
        if ~exist('option','var'), return; end
        
        switch option
            case 'view'
                viewTwix(h);
            case 'read'
                setWaitbar(h);
                readTwix(h);
            case 'map'
                setWaitbar(h);
                mapTwix(h)
        end

    end % end of constructor
        
    function open(h, fname, SelectedPart)
    %% open(fname) -- constructor of getwix open TWIX file and protocol
    %
    % Allows to read several twix files while keeping the 'pub' property
    %
    % See also: setBuffer, setWaitbar, fetchReadout, close
    
        if h.fid ~= -1
            disp( [ mfilename('class') ': file not closed: ' h.File.name ]);
            return
        end

        % check filename
        if ~exist(fname,'file')
            disp( [ mfilename('class') ': file not found: ' fname ]);
            return
        end
            
        % check file type
        d = dir(fname);
        if length(d) > 1 || d(1).isdir
            disp( [ mfilename('class') ': filename is directory: ' fname ]);
            return
        end
        
        % header and protocol
        [h.fid, errmsg] = fopen(fname,'r');
        if h.fid < 0
            h.fid      = -1;
            disp( [ mfilename('class') ': ' errmsg ] );
            return
        end
        
        h.File.name     = fname;
        [h.File.fullname, h.File.permission, h.File.machinefmt] = fopen(h.fid);
        [h.File.path, h.File.basename, e]=fileparts(h.File.fullname);
        h.File.syngo  	= h.syngoVB; % default structure for VB15, VB17, VB19
        h.File.name     = [ h.File.basename e ];
        h.File.sizeTwix = d.bytes;
        h.File.size     = d.bytes;
        h.File.FieldNames  = {};
        
        % VB and VD structure
        h.File.oskip    = fread(h.fid, 1, 'ulong');
        h.File.pos      = 4;
        
        % guess VD format
        h.File.selectedPart = 1;
        
        if h.File.oskip == 0
            h.File.syngo    = h.syngoVD; % structure for VD
            n = fread(h.fid, 1, 'ulong');
            for k=1:n
                kk = fread(h.fid, 6, 'ulong');
                h.File.twixpart(k).oskip = kk(3);
                h.File.twixpart(k).size = kk(5);
                st = fread(h.fid,16*8,'uint8=>char')';
                h.File.twixpart(k).sequence = deblank(st(65:end));
                if n>1
                    fprintf(1,'part %d: %s (oskip=%d size=%d)\n', k, ...
                        h.File.twixpart(k).sequence, ...
                        h.File.twixpart(k).oskip, ...
                        h.File.twixpart(k).size);
                end
            end
            if SelectedPart == -1, SelectedPart = n; end
            h.File.selectedPart = SelectedPart;
            fseek(h.fid, h.File.twixpart(h.File.selectedPart).oskip, -1);
            h.File.oskip = fread(h.fid, 1, 'ulong');
            h.File.oskip = h.File.oskip + h.File.twixpart(k).oskip;
            h.File.pos   = h.File.pos + h.File.twixpart(k).oskip;
            h.File.size  = h.File.twixpart(k).oskip + h.File.twixpart(k).size;
        end
        
        if h.File.oskip > h.File.size
            fclose(h.fid);
            h.fid = -1;
            disp([ mfilename('class') ': not twix format: ' fname ]);
            return
        end
        
        % headers
        h.XProt.hdr = fread(h.fid, h.File.oskip - h.File.pos, 'uint8=>char')';
        h.File.pos = h.File.oskip;
        
       if strcmp(e,'.out')
            
            % VA35 structure
            h.File.syngo  	= h.syngoVA;
            fname = fullfile(h.File.path, [h.File.basename '.asc']);
            if exist(fname, 'file')
                fp = fopen(fname, 'r');
                hdr = fread(fp, inf, 'uint8=>char')';
                fclose(fp);
                fetchMrProt(h, hdr);
                
                % add WIP MemBlock from usertwix
                try
                    fhdl = h.MrProt.tSequenceFileName;
                    fhdl = fhdl(max([ strfind(fhdl,'\') strfind(fhdl,'/') ])+1:end);
                    fhdl = str2func( [ 'usertwix' fhdl ] );
                    h.File.wip = fhdl(true);
                catch
                    h.File.wip = [];
                end
            end
            
            fname = fullfile(h.File.path, 'EVA95.tmp');
            if exist(fname, 'file')
                fp = fopen(fname, 'r');
                hdr = fread(fp, inf, 'uint8=>char')';
                fclose(fp);
                fetchMrXML(h, hdr)
                dum.EVA = h.XProt;
                h.XProt.EVA = dum.EVA;
            end
            
       else
        
           if fetchMrProt(h, h.XProt.hdr)
               
               % add WIP MemBlock from usertwix (if available)
               try
                   fhdl = h.MrProt.tSequenceFileName;
                   fhdl = fhdl(max([ strfind(fhdl,'\') strfind(fhdl,'/') ])+1:end);
                   fhdl = str2func( [ 'usertwix' fhdl ] );
                   h.File.wip = fhdl(true);
               catch
                   h.File.wip = [];
               end
           else
               disp( [ mfilename('class') ': error while reading MrProt: ' fname ] );
           end
           
       end

        
        % read the first MDH to determine some initial parameters
        DMAlen1                 = fread(h.fid, 1, 'ushort');
        DMAlen2                 = fread(h.fid, 2, 'uint8');
        DMAlength               = DMAlen1 + DMAlen2(1) * 65536;
        
        if feof(h.fid) || DMAlength == 0   || DMAlength > h.File.size
            disp( [ mfilename('class') ': first DMA not found in: ' fname ] );
            close(h);
            return
        end
        
        h.File.MeasUID          = fread(h.fid, 1, 'long');
        h.BufferMxByte          = setBuffer(h);
        h.BufferMxByte          = ceil(h.BufferMxByte / DMAlength) * DMAlength ;
        
        fseek(h.fid, 4, 'cof');
        h.File.TimeStampInHour    = 2.5 * double(fread(h.fid, 1, 'ulong')) / 1000 / 60 / 60;
        h.File.PMUTimeStampInHour = 2.5 * double(fread(h.fid, 1, 'ulong')) / 1000 / 60 / 60;
        
        % rewind to oskip position of twix
        fseek(h.fid, h.File.oskip, -1);
        
    end % end of open

    function map(h)
        % maps twix file in h.pub.map array
        % description of dimension is stored in h.pub.mapStr
        h.close();
        h.open(h.File.fullname, h.File.selectedPart);
        mapTwix(h);
    end
        
    function getXProt(h)
        if isfield(h.XProt, 'hdr')
            fetchMrXML(h, h.XProt.hdr);
        end
    end
    
    function close(h, CancelRequested)
    %% close twix file (without deleteing the pub properties)
    %
    % See also: open
    
        if h.fid > -1
            if nargin>1 && CancelRequested, h.fid = -2; end
            if ~isempty(h.WaitBar), delete(h.WaitBar); end
            fclose(h.fid);
            h.Buffer       = [];
            h.fid          = -1;
            h.BufferedDMA  = false;
            h.BufferPos    = 0;
            h.BufferSiz    = 0;
            h.BufferMxByte = 0;
            h.WithWaitBar  = -1;
            h.WaitBar      = [];
        end
    end % of close

    function status = ok(h)
    %% ok() -- return true when current file is opened
    %
    % See also: fetchReadout, close, open
        status = ( (h.fid > 0) &&  strcmp(h.File.fullname, fopen(h.fid)) );
    end % of ok
    
    function bufferSize = setBuffer(h, Value)
    %% setBuffer( false or SizeInMB ) -- default buffer size 8 MBytes
    %
    % DMA length (MrServers\MrHardwareControl\PciDsp\RX\RXInstr.h)
    %
    %  For each of PCI_RX board #0..3 of MRIR #0:
    %  The following information is encoded:
    %   1. DMA-Length
    %   2. Pack Bit
    %   3. Boolean flags if a board is used
    %     
    %   Bit
    %    31     28 27  26 25 24                      0
    %   +---------+------+--+-------------------------+
    %   | <enable>|unused|  |   <DMA length>          |
    %   +---------+------+--+-------------------------+
    %    \---+---/        \/ \--------+--------------/
    %        |            |           |
    %        |            |  DMA length in units of bytes:
    %        |            |  = Number of used channels of this PCI_RX board
    %        |            |    * (MDH size + echo length * 2 * 4)
    %        |            |                                |   |
    %        |            |                             Re/Im  4 byte/float
    %        |            +-- Pack Bit:
    %        |                0: 1st readout of DMA
    %        |                1: 2nd por following readout of DMA
    %        |
    %        +--------------- PCI_RX board enable flags:
    %                         Bit 28:  PCI_RX 3 enable
    %                         Bit 29:  PCI_RX 2 enable
    %                         Bit 30:  PCI_RX 1 enable
    %                         Bit 31:  PCI_RX 0 enable
    %
    % See also: fetchReadout, open, close
    
        if h.BufferSiz ~= 0
            disp( [ mfilename('class') ': attempt to change buffered DMA while reading: ' h.File.name ]);
            
        elseif nargin == 1
            h.BufferedDMA  = true;
            h.BufferMxByte = 8 * 1024 * 1024; % default 8 MB
            
        elseif Value > 0
            h.BufferedDMA  = true;
            h.BufferMxByte = round(Value * 1024 * 1024);
               
        else
            h.BufferedDMA  = false;
            h.BufferMxByte = 0;
        end
        
        if nargout == 1, bufferSize = h.BufferMxByte; end
    end
    
    
    function nCoils = fetchReadout(h)
    %% fetchReadout() -- get Data Measurement Array
    %  syngo VB15 & VE11 Measurement Data Header (Mdh) + channels * ADC
	%
    %% VB15 structure
    % All array channels for the current readout are read, including the  measurement data header (128 bytes long)
    % and ADC (complex). These are stored in Mdh(CoilNumber).
    %
    %               DMAlength: uint32 ulFlagsAndDMALength           // bit  0..24: DMA length bytes
    %                DMAflags:                                      // bit 25: pack bit, bit 26..31: pci_rx enable flags
    %                 MeasUID: int32  lMeasUID                      // measurement user ID
    %             ScanCounter: uint32 ulScanCounter                 // scan counter [1...]
    %               TimeStamp: uint32 ulTimeStamp                   // time stamp [2.5 ms ticks since 00:00]
    %           TimeStampInMs:
    %            PMUTimeStamp: uint32 ulPMUTimeStamp                // PMU time stamp [2.5 ms ticks since last trigger]
    %        PMUTimeStampInMs:
    %            EvalInfoMask: uint32 aulEvalInfoMask[MDH_NUMBEROFEVALINFOMASK] // evaluation info mask field
    %         EvalInfoMaskBit:
    %           EvalInfoMask2: uint32 aulEvalInfoMask[MDH_NUMBEROFEVALINFOMASK] // 2nd evaluation info mask field
    %           SamplesInScan: uint16 ushSamplesInScan              // # of samples acquired in scan
    %            UsedChannels: uint16 ushUsedChannels               // # of channels used in scan
    %                    Clin: uint16 ushLine                       // line index
    %                    Cacq: uint16 ushAcquisition                // acquisition index
    %                    Cslc: uint16 ushSlice                      // slice index
    %                    Cpar: uint16 ushPartition                  // partition index
    %                    Ceco: uint16 ushEcho                       // echo index
    %                    Cphs: uint16 ushPhase                      // phase index
    %                    Crep: uint16 ushRepetition                 // measurement repeat index
    %                    Cset: uint16 ushSet                        // set index
    %                    Cseg: uint16 ushSeg                        // segment index  (for TSE)
    %                    Cida: uint16 ushIda                        // IceDimension a index
    %                    Cidb: uint16 ushIdb                        // IceDimension b index
    %                    Cidc: uint16 ushIdc                        // IceDimension c index
    %                    Cidd: uint16 ushIdd                        // IceDimension d index
    %                    Cide: uint16 ushIde                        // IceDimension e index
    %           CutOffDataPre: uint16 sCutOffData.Pre               // cut-off values
    %          CutOffDataPost: uint16 sCutOffData.Post
    %      KSpaceCentreColumn: uint16 ushKSpaceCentreColumn         // centre of echo
    %              CoilSelect: uint16 ushCoilSelect                 // bit 0..3: CoilSelect
    %        ReadOutOffcentre: float  fReadOutOffcentre             // ReadOut offcenter value
    %         TimeSinceLastRF: uint32 ulTimeSinceLastRF             // Sequence time stamp since last RF pulse
    %      KSpaceCentreLineNo: uint16 ushKSpaceCentreLineNo         // number of K-space centre line
    % KSpaceCentrePartitionNo: uint16 ushKSpaceCentrePartitionNo    // number of K-space centre partition
    %          IceProgramPara: uint16 aushIceProgramPara[MDH_NUMBEROFICEPROGRAMPARA] // free parameter for IceProgram
    %                FreePara: uint16 aushFreePara[MDH_FREEHDRPARA] // free parameter
    %                      SD: sSliceData                           // 28 bytes
    %               ChannelId: uint16 ushChannelId                  // channel Id must be the last parameter
    %              PTABPosNeg: uint16 ushPTABPosNeg                 // negative, absolute PTAB position in [0.1 mm]
    %
    %                     adc: double
    %
    %% VE11 structure:
	% // Source: /n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h
	% // Definition of measurement data header -- VE11 -- total length: 6 x 32 Byte (192 Byte)
	% #define MDH_NUMBEROFEVALINFOMASK   2
	% #define MDH_NUMBEROFICEPROGRAMPARA 24
	% #define MDH_RESERVEDHDRPARA  (4)
	% 
	% // Definition of time stamp tick interval/frequency                         */
	% #define RXU_TIMER_INTERVAL  (2500000)     /* data header timer interval [ns]*/
	%
	% PACKED_MEMBER( uint32_t,     ulFlagsAndDMALength           );                 ///<  0: ( 4) bit  0..24: DMA length [bytes]
	%																			    ///<          bit     25: pack bit
	%																			    ///<          bit 26..31: pci_rx enable flags
	% PACKED_MEMBER( int32_t,      lMeasUID                      );                 ///<  4: ( 4) measurement user ID
	% PACKED_MEMBER( uint32_t,     ulScanCounter                 );                 ///<  8: ( 4) scan counter [1...]
	% PACKED_MEMBER( uint32_t,     ulTimeStamp                   );                 ///< 12: ( 4) time stamp [2.5 ms ticks since 00:00]
	% PACKED_MEMBER( uint32_t,     ulPMUTimeStamp                );                 ///< 16: ( 4) PMU time stamp [2.5 ms ticks since last trigger]
	% PACKED_MEMBER( uint16_t,     ushSystemType                 );                 ///< 20: ( 2) System type (todo: values?? ####)
	% PACKED_MEMBER( uint16_t,     ulPTABPosDelay                );                 ///< 22: ( 2z) PTAb delay ??? TODO: How do we handle this ####
	% PACKED_MEMBER( int32_t,	     lPTABPosX                     );                 ///< 24: ( 4) absolute PTAB position in [µm]
	% PACKED_MEMBER( int32_t,	     lPTABPosY                     );                 ///< 28: ( 4) absolute PTAB position in [µm]
	% PACKED_MEMBER( int32_t,	     lPTABPosZ                     );                 ///< 32: ( 4) absolute PTAB position in [µm]
	% PACKED_MEMBER( uint32_t,	   ulReserved1                   );                 ///< 36: ( 4) reserved for future hardware signals
	% PACKED_MEMBER( uint32_t,     aulEvalInfoMask[MDH_NUMBEROFEVALINFOMASK]);      ///< 40: ( 8) evaluation info mask field
	% PACKED_MEMBER( uint16_t,     ushSamplesInScan              );                 ///< 48: ( 2) # of samples acquired in scan
	% PACKED_MEMBER( uint16_t,     ushUsedChannels               );                 ///< 50: ( 2) # of channels used in scan
	% PACKED_STRUCT( sLoopCounter, sLC                           );                 ///< 52: (28) loop counter structure
	%   PACKED_MEMBER( uint16_t,  ushLine         ); /**< line index                   */
	%   PACKED_MEMBER( uint16_t,  ushAcquisition  ); /**< acquisition index            */
	%   PACKED_MEMBER( uint16_t,  ushSlice        ); /**< slice index                  */
	%   PACKED_MEMBER( uint16_t,  ushPartition    ); /**< partition index              */
	%   PACKED_MEMBER( uint16_t,  ushEcho         ); /**< echo index                   */
	%   PACKED_MEMBER( uint16_t,  ushPhase        ); /**< phase index                  */
	%   PACKED_MEMBER( uint16_t,  ushRepetition   ); /**< measurement repeat index     */
	%   PACKED_MEMBER( uint16_t,  ushSet          ); /**< set index                    */
	%   PACKED_MEMBER( uint16_t,  ushSeg          ); /**< segment index  (for TSE)     */
	%   PACKED_MEMBER( uint16_t,  ushIda          ); /**< IceDimension a index         */
	%   PACKED_MEMBER( uint16_t,  ushIdb          ); /**< IceDimension b index         */
	%   PACKED_MEMBER( uint16_t,  ushIdc          ); /**< IceDimension c index         */
	%   PACKED_MEMBER( uint16_t,  ushIdd          ); /**< IceDimension d index         */
	%   PACKED_MEMBER( uint16_t,  ushIde          ); /**< IceDimension e index         */
	% PACKED_STRUCT( sCutOffData,  sCutOff                       );                 ///< 80: ( 4) cut-off values
    %   PACKED_MEMBER( uint16_t,  ushPre          );    /**< write ushPre zeros at line start */
    %   PACKED_MEMBER( uint16_t,  ushPost         );    /**< write ushPost zeros at line end  */
	% PACKED_MEMBER( uint16_t,     ushKSpaceCentreColumn         );                 ///< 84: ( 2) centre of echo
	% PACKED_MEMBER( uint16_t,     ushCoilSelect                 );                 ///< 86: ( 2) Bit 0..3: CoilSelect
	% PACKED_MEMBER( float,        fReadOutOffcentre             );                 ///< 88: ( 4) ReadOut offcenter value
	% PACKED_MEMBER( uint32_t,     ulTimeSinceLastRF             );                 ///< 92: ( 4) Sequence time stamp since last RF pulse
	% PACKED_MEMBER( uint16_t,     ushKSpaceCentreLineNo         );                 ///< 96: ( 2) number of K-space centre line
	% PACKED_MEMBER( uint16_t,     ushKSpaceCentrePartitionNo    );                 ///< 98: ( 2) number of K-space centre partition
	% PACKED_STRUCT( sSliceData,   sSD                           );                 ///< 100:(28) Slice Data
    %   PACKED_STRUCT( sVector,         sSlicePosVec     ); /**< slice position vector        */
    %     PACKED_MEMBER( float,  flSag          );
    %     PACKED_MEMBER( float,  flCor          );
    %     PACKED_MEMBER( float,  flTra          );
    %   PACKED_MEMBER( float,           aflQuaternion[4] ); /**< rotation matrix as quaternion*/
	% PACKED_MEMBER( uint16_t,     aushIceProgramPara[MDH_NUMBEROFICEPROGRAMPARA] );///< 128:(48) free parameter for IceProgram
	% PACKED_MEMBER( uint16_t,     aushReservedPara[MDH_RESERVEDHDRPARA] );         ///< 176:( 8) unused parameter (padding to next 192byte alignment )
	%																			    ///<          NOTE: These parameters MUST NOT be used by any application (for future use)
	% PACKED_MEMBER( uint16_t,     ushApplicationCounter         );                 ///< 184 ( 2)
	% PACKED_MEMBER( uint16_t,     ushApplicationMask            );                 ///< 186 ( 2)
	% PACKED_MEMBER( uint32_t,     ulCRC                         );                 ///< 188:( 4) CRC 32 checksum
    %
	%% channel data header
	%  PACKED_MEMBER( uint32_t,     ulTypeAndChannelLength        );    ///< 0: (4) bit  0.. 7: type (0x02 => ChannelHeader)
	%                                                                   ///<        bit  8..31: channel length (header+data) in byte
	%                                                                   ///<        type   := ulTypeAndChannelLength & 0x000000FF
	%                                                                   ///<        length := ulTypeAndChannelLength >> 8
	%  PACKED_MEMBER( int32_t,      lMeasUID                      );    ///< 4: (4) measurement user ID
	%  PACKED_MEMBER( uint32_t,     ulScanCounter                 );    ///< 8: (4) scan counter [1...]
	%  PACKED_MEMBER( uint32_t,     ulReserved1                   );    ///< 12:(4) reserved
	%  PACKED_MEMBER( uint32_t,     ulSequenceTime                );    ///< 16:(4) Sequence readout starting time bit 31..9 time in [10us]
	%                                                                   ///<                                       bit  8..0 time in [25ns]
	%  PACKED_MEMBER( uint32_t,     ulUnused2                     );    ///< 20:(4) unused
	%  PACKED_MEMBER( uint16_t,     ulChannelId                   );    ///< 24:(4) unused
	%  PACKED_MEMBER( uint16_t,     ulUnused3                     );    ///< 26:(2) unused
	%  PACKED_MEMBER( uint32_t,     ulCRC                         );    ///< 28:(4) CRC32 checksum of channel header
    %
    % See also: setBuffer, infoMDH
    
        nCoils = 0;
        h.Cdh = []; % reset channels
        
        if h.fid < 0
            disp( [ mfilename('class') ': no current twix file' ] );
            if h.fid == -2, error('cancelling ...'); end
            return;
        end
        
        if h.BufferedDMA
            nCoils = getReadoutByBuffer(h);
        else
            nCoils = getReadoutByRead(h);
        end
        
        if MDH_SYNCDATA(h)
            nCoils = 1;
            fprintf(1,'SYNCDATA detected');
        end
        
        if MDH_ACQEND(h)
            nCoils = 0;
        end
        
        updateWaitbar(h);
        
    end % end of fetchReadout
    
    %% Loop counters (Matlab convention)

    function val = lin(h)
    % lin loop counter (h.Mdh.Loop.Clin + 1)
        val = 1 + h.Mdh.Loop.Clin;
    end
    function val = acq(h)
    % acq loop counter (h.Mdh.Loop.Cacq + 1)
        val = 1 + h.Mdh.Loop.Cacq;
    end
    function val = slc(h)
    % slc loop counter (h.Mdh.Loop.Cslc + 1)
        val = 1 + h.Mdh.Loop.Cslc;
    end
    function val = par(h)
    % par loop counter (h.Mdh.Loop.Cpar + 1)
        val = 1 + h.Mdh.Loop.Cpar;
    end
    function val = eco(h)
    % eco loop counter (h.Mdh.Loop.Ceco + 1)
        val = 1 + h.Mdh.Loop.Ceco;
    end
    function val = phs(h)
    % phs loop counter (h.Mdh.Loop.Cphs + 1)
        val = 1 + h.Mdh.Loop.Cphs;
    end
    function val = rep(h)
    % rep loop counter (h.Mdh.Loop.Crep + 1)
        val = 1 + h.Mdh.Loop.Crep;
    end
    function val = set(h)
    % set loop counter (h.Mdh.Loop.Cset + 1)
        val = 1 + h.Mdh.Loop.Cset;
    end
    function val = seg(h)
    % seg loop counter (h.Mdh.Loop.Cseg + 1)
        val = 1 + h.Mdh.Loop.Cseg;
    end
    function val = ida(h)
    % ida loop counter (h.Mdh.Loop.Cida + 1)
        val = 1 + h.Mdh.Loop.Cida;
    end
    function val = idb(h)
    % idb loop counter (h.Mdh.Loop.Cidb + 1)
        val = 1 + h.Mdh.Loop.Cidb;
    end
    function val = idc(h)
    % idc loop counter (h.Mdh.Loop.Cidc + 1)
        val = 1 + h.Mdh.Loop.Cidc;
    end
    function val = idd(h)
    % idd loop counter (h.Mdh.Loop.Cidd + 1)
        val = 1 + h.Mdh.Loop.Cidd;
    end
    function val = ide(h)
    % ide loop counter (h.Mdh.Loop.Cide + 1)
        val = 1 + h.Mdh.Loop.Cide;
    end

    %% EvalInfoMask (Matlab bitget convention)
    
    function Flag = MDH_ACQEND(h)
    % MDH_ACQEND (1) -- last DMA signals end of acquisition
        Flag = bitget(h.Mdh.EvalInfoMask, 1);
    end
    function Flag = MDH_RTFEEDBACK(h)
    % MDH_RTFEEDBACK (2) -- priority flag for real time feedback scan
        Flag = bitget(h.Mdh.EvalInfoMask, 2);
    end
    function Flag = MDH_HPFEEDBACK(h)
    % MDH_HPFEEDBACK (3) -- high priority real time
        Flag = bitget(h.Mdh.EvalInfoMask, 3);
    end
    function Flag = MDH_ONLINE(h)
    % MDH_ONLINE (4) -- processed by ICE
        Flag = bitget(h.Mdh.EvalInfoMask, 4);
    end
    function Flag = MDH_OFFLINE(h)
    % MDH_OFFLINE (5) -- ignored by ICE
        Flag = bitget(h.Mdh.EvalInfoMask, 5);
    end
    function Flag = MDH_SYNCDATA(h)
    % MDH_SYNCDATA (6) -- readout contains synchroneous data
        Flag = bitget(h.Mdh.EvalInfoMask, 6);
    end
    function Flag = MDH_LASTSCANINCONCAT(h)
    % MDH_LASTSCANINCONCAT (9) -- Flag for last scan in concatenation
        Flag = bitget(h.Mdh.EvalInfoMask, 9);
    end
    function Flag = MDH_RAWDATACORRECTION(h)
    % MDH_RAWDATACORRECTION (11) -- Correct rawadata with rawdata correction factor
        Flag = bitget(h.Mdh.EvalInfoMask, 11);
    end
    function Flag = MDH_LASTSCANINMEAS(h)
    % MDH_LASTSCANINMEAS (12) -- Flag for last scan in measurement
        Flag = bitget(h.Mdh.EvalInfoMask, 12);
    end
    function Flag = MDH_SCANSCALEFACTOR(h)
    % MDH_SCANSCALEFACTOR (13) -- Flag for scan specific additional scale factor
        Flag = bitget(h.Mdh.EvalInfoMask, 13);
    end
    function Flag = MDH_2NDHADAMARPULSE(h)
    % MDH_2NDHADAMARPULSE (14) -- 2nd RF exitation of HADAMAR
        Flag = bitget(h.Mdh.EvalInfoMask, 14);
    end
    function Flag = MDH_REFPHASESTABSCAN(h)
    % MDH_REFPHASESTABSCAN (15) -- reference phase stabilization scan
        Flag = bitget(h.Mdh.EvalInfoMask, 15);
    end
    function Flag = MDH_PHASESTABSCAN(h)
    % MDH_PHASESTABSCAN (16) - phase stabilization scan
        Flag = bitget(h.Mdh.EvalInfoMask, 16);
    end
    function Flag = MDH_D3FFT(h)
    % MDH_D3FFT (17) -- execute 3D FFT
        Flag = bitget(h.Mdh.EvalInfoMask, 17);
    end
    function Flag = MDH_SIGNREV(h)
    % MDH_SIGNREV (18) -- sign reversal
        Flag = bitget(h.Mdh.EvalInfoMask, 18);
    end
    function Flag = MDH_PHASEFFT(h)
    % MDH_PHASEFFT (19) -- execute phase fft
        Flag = bitget(h.Mdh.EvalInfoMask, 19);
    end
    function Flag = MDH_SWAPPED(h)
    % MDH_SWAPPED (20) -- swapped phase/readout direction
        Flag = bitget(h.Mdh.EvalInfoMask, 20);
    end
    function Flag = MDH_POSTSHAREDLINE(h)
    % MDH_POSTSHAREDLINE (21) -- shared line
        Flag = bitget(h.Mdh.EvalInfoMask, 21);
    end
    function Flag = MDH_PHASCOR(h)
    % MDH_PHASCOR (22) -- phase correction data
        Flag = bitget(h.Mdh.EvalInfoMask, 22);
    end
    function Flag = MDH_PATREFSCAN(h)
    % MDH_PATREFSCAN (23) -- additonal scan for PAT reference line/partition
        Flag = bitget(h.Mdh.EvalInfoMask, 23);
    end
    function Flag = MDH_PATREFANDIMASCAN(h)
    % MDH_PATREFANDIMASCAN (24) -- additonal scan for PAT reference line/partition that is also used as image scan
        Flag = bitget(h.Mdh.EvalInfoMask, 24);
    end
    function Flag = MDH_REFLECT(h)
    % MDH_REFLECT (25) -- reflect line
        Flag = bitget(h.Mdh.EvalInfoMask, 25);
    end
    function Flag = MDH_NOISEADJSCAN(h)
    % MDH_NOISEADJSCAN (26) -- noise adjust scan
        Flag = bitget(h.Mdh.EvalInfoMask, 26);
    end
    function Flag = MDH_SHARENOW(h)
    % MDH_SHARENOW (27) -- all lines are acquired from the actual and previous e.g. phases
        Flag = bitget(h.Mdh.EvalInfoMask, 27);
    end
    function Flag = MDH_LASTMEASUREDLINE(h)
    % MDH_LASTMEASUREDLINE (28) -- indicates that current line is last measured line of all succeeding e.g. phases
        Flag = bitget(h.Mdh.EvalInfoMask, 28);
    end
    function Flag = MDH_FIRSTSCANINSLICE(h)
    % MDH_FIRSTSCANINSLICE (29) -- indicates first scan in slice (needed for time stamps)
        Flag = bitget(h.Mdh.EvalInfoMask, 29);
    end
    function Flag = MDH_LASTSCANINSLICE(h)
    % MDH_LASTSCANINSLICE (30) -- indicates last scan in slice (needed for time stamps)
        Flag = bitget(h.Mdh.EvalInfoMask, 30);
    end
    function Flag = MDH_TREFFECTIVEBEGIN(h)
    % MDH_TREFFECTIVEBEGIN (31) -- indicates begin time stamp for TReff (triggered measurement)
        Flag = bitget(h.Mdh.EvalInfoMask, 31);
    end
    function Flag = MDH_TREFFECTIVEEND(h)
    % MDH_TREFFECTIVEEND (32) -- indicates end time stamp for TReff (triggered measurement)
        Flag = bitget(h.Mdh.EvalInfoMask, 32);
    end
   
    %% EvalInfoMask2 (Matlab bitget convention)

    function Flag = MDH_MDS_REF_POSITION(h)
    % MDH_MDS_REF_POSITION (32 + 1) -- indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
        Flag = bitget(h.Mdh.EvalInfoMask2, 1);
    end
    function Flag = MDH_SLC_AVERAGED(h)
    % MDH_SLC_AVERAGED (32 + 2) -- indicates avveraged slice for slice partial averaging scheme
        Flag = bitget(h.Mdh.EvalInfoMask2, 2);
    end
    function Flag = MDH_TAGFLAG1(h)
    % MDH_TAGFLAG1 (32 + 3) -- adjust scan
        Flag = bitget(h.Mdh.EvalInfoMask2, 3);
    end
    function Flag = MDH_CT_NORMALIZE(h)
    % MDH_CT_NORMALIZE (32 + 4) -- Marks scans used to calculate correction maps for TimCT-Prescan normalize
        Flag = bitget(h.Mdh.EvalInfoMask2, 4);
    end
    function Flag = MDH_SCAN_FIRST(h)
    % MDH_SCAN_FIRST (32 + 5) -- Marks the first scan of a particular map
        Flag = bitget(h.Mdh.EvalInfoMask2, 5);
    end
    function Flag = MDH_SCAN_LAST(h)
    % MDH_SCAN_LAST (32 + 6) -- Marks the last scan of a particular map
        Flag = bitget(h.Mdh.EvalInfoMask2, 6);
    end
    function Flag = MDH_FIRST_SCAN_IN_BLADE(h)
    % MDH_FIRST_SCAN_IN_BLADE (32 + 9) -- Marks the first line of a blade
        Flag = bitget(h.Mdh.EvalInfoMask2, 9);
    end
    function Flag = MDH_LAST_SCAN_IN_BLADE(h)
    % MDH_LAST_SCAN_IN_BLADE (32 + 10) -- Marks the last line of a blade
        Flag = bitget(h.Mdh.EvalInfoMask2, 10);
    end
    function Flag = MDH_LAST_BLADE_IN_TR(h)
    % MDH_LAST_BLADE_IN_TR (32 + 11) -- Set for all lines of the last BLADE in each TR interval
        Flag = bitget(h.Mdh.EvalInfoMask2, 11);
    end    
    function Flag = MDH_PACE(h)
    % MDH_PACE (32 + 13) -- Distinguishes PACE scans from non PACE scans
        Flag = bitget(h.Mdh.EvalInfoMask2, 13);
    end
    function Flag = MDH_RETRO_LASTPHASE(h)
    % MDH_RETRO_LASTPHASE (32 + 14) -- Marks the last phase in a heartbeat
        Flag = bitget(h.Mdh.EvalInfoMask2, 14);
    end
    function Flag = MDH_RETRO_ENDOFMEAS(h)
    % MDH_RETRO_ENDOFMEAS (32 + 15) -- Marks an ADC at the end of the measurement
        Flag = bitget(h.Mdh.EvalInfoMask2, 15);
    end
    function Flag = MDH_RETRO_REPEATTHISHEARTBEAT(h)
    % MDH_RETRO_REPEATTHISHEARTBEAT (32 + 16) -- Repeat the current heartbeat when this bit is found
        Flag = bitget(h.Mdh.EvalInfoMask2, 16);
    end
    function Flag = MDH_RETRO_REPEATPREVHEARTBEAT(h)
    % MDH_RETRO_REPEATPREVHEARTBEAT (32 + 17) -- Repeat the previous heartbeat when this bit is found
        Flag = bitget(h.Mdh.EvalInfoMask2, 17);
    end
    function Flag = MDH_RETRO_ABORTSCANNOW(h)
    % MDH_RETRO_ABORTSCANNOW (32 + 18) -- Just abort everything
        Flag = bitget(h.Mdh.EvalInfoMask2, 18);
    end
    function Flag = MDH_RETRO_LASTHEARTBEAT(h)
    % MDH_RETRO_LASTHEARTBEAT (32 + 19) -- This adc is from the last heartbeat (a dummy)
        Flag = bitget(h.Mdh.EvalInfoMask2, 19);
    end
    function Flag = MDH_RETRO_DUMMYSCAN(h)
    % MDH_RETRO_DUMMYSCAN (32 + 20) -- This adc is just a dummy scan, throw it away
        Flag = bitget(h.Mdh.EvalInfoMask2, 20);
    end
    function Flag = MDH_RETRO_ARRDETDISABLED(h)
    % MDH_RETRO_ARRDETDISABLED (32 + 21) -- Disable all arrhythmia detection when this bit is found
        Flag = bitget(h.Mdh.EvalInfoMask2, 21);
    end
    function Flag = MDH_B1_CONTROLLOOP(h)
    % MDH_B1_CONTROLLOOP (32 + 22) -- Marks the readout as to be used for B1 Control Loop
        Flag = bitget(h.Mdh.EvalInfoMask2, 22);
    end
    function Flag = MDH_SKIP_ONLINE_PHASCOR(h)
    % MDH_SKIP_ONLINE_PHASCOR (32 + 23) -- Marks scans not to be online phase corrected, even if online phase correction is switched on
        Flag = bitget(h.Mdh.EvalInfoMask2, 23);
    end
    function Flag = MDH_SKIP_REGRIDDING(h)
    % MDH_SKIP_REGRIDDING (32 + 24) -- Marks scans not to be regridded, even if regridding is switched on
        Flag = bitget(h.Mdh.EvalInfoMask2, 24);
    end    

    function MDH = infoMDH(h)
    %% strInfoMask(MdhNumber) - evaluates Mdh Info Mask as string
    %
    % EvalInfoMask
    % MDH_ACQEND            (0UL)   last scan
    % MDH_RTFEEDBACK        (1UL)   Realtime feedback scan
    % MDH_HPFEEDBACK        (2UL)   High perfomance feedback scan
    % MDH_ONLINE            (3UL)   processing should be done online
    % MDH_OFFLINE           (4UL)   processing should be done offline
    % MDH_SYNCDATA          (5UL)   readout contains synchroneous data
    % MDH_LASTSCANINCONCAT  (8UL)   Flag for last scan in concatination
    % MDH_RAWDATACORRECTION (10UL)   Correct the rawadata with the rawdata correction factor
    % MDH_LASTSCANINMEAS    (11UL)   Flag for last scan in measurement
    % MDH_SCANSCALEFACTOR   (12UL)   Flag for scan specific additional scale factor
    % MDH_2NDHADAMARPULSE   (13UL)   2nd RF exitation of HADAMAR
    % MDH_REFPHASESTABSCAN  (14UL)   reference phase stabilization scan
    % MDH_PHASESTABSCAN     (15UL)   phase stabilization scan
    % MDH_D3FFT             (16UL)   execute 3D FFT
    % MDH_SIGNREV           (17UL)   sign reversal
	% MDH_PHASEFFT          (18UL)   execute phase fft
    % MDH_SWAPPED           (19UL)   swapped phase/readout direction
    % MDH_POSTSHAREDLINE    (20UL)   shared line
    % MDH_PHASCOR           (21UL)   phase correction data
    % MDH_PATREFSCAN        (22UL)   additonal scan for PAT reference line/partition
    % MDH_PATREFANDIMASCAN  (23UL)   additonal scan for PAT reference line/partition that is also used as image scan
    % MDH_REFLECT           (24UL)   reflect line
    % MDH_NOISEADJSCAN      (25UL)   noise adjust scan
    % MDH_SHARENOW          (26UL)   all lines are acquired from the actual and previous e.g. phases
    % MDH_LASTMEASUREDLINE  (27UL)   indicates that the current line is the last measured line of all succeeding e.g. phases
    % MDH_FIRSTSCANINSLICE  (28UL)   indicates first scan in slice (needed for time stamps)
    % MDH_LASTSCANINSLICE   (29UL)   indicates  last scan in slice (needed for time stamps)
    % MDH_TREFFECTIVEBEGIN  (30UL)   indicates the begin time stamp for TReff (triggered measurement)
    % MDH_TREFFECTIVEEND    (31UL)   indicates the   end time stamp for TReff (triggered measurement)
    %
    % EvalInfoMask2
    % MDH_MDS_REF_POSITION          (32UL)  indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
    % MDH_SLC_AVERAGED              (33UL)  indicates avveraged slice for slice partial averaging scheme
    % MDH_TAGFLAG1                  (34UL)  adjust scan
    % MDH_CT_NORMALIZE              (35UL)  Marks scans used to calculate correction maps for TimCT-Prescan normalize
    % MDH_SCAN_FIRST                (36UL)  Marks the first scan of a particular map
    % MDH_SCAN_LAST                 (37UL)  Marks the last scan of a particular map
    % MDH_FIRST_SCAN_IN_BLADE       (40UL)  Marks the first line of a blade
    % MDH_LAST_SCAN_IN_BLADE        (41UL)  Marks the last line of a blade
    % MDH_LAST_BLADE_IN_TR          (42UL)  Set for all lines of the last BLADE in each TR interval
    % MDH_PACE                      (44UL)  Distinguishes PACE scans from non PACE scans
    % MDH_RETRO_LASTPHASE           (45UL)  Marks the last phase in a heartbeat
    % MDH_RETRO_ENDOFMEAS           (46UL)  Marks an ADC at the end of the measurement
    % MDH_RETRO_REPEATTHISHEARTBEAT (47UL)  Repeat the current heartbeat when this bit is found
    % MDH_RETRO_REPEATPREVHEARTBEAT (48UL)  Repeat the previous heartbeat when this bit is found
    % MDH_RETRO_ABORTSCANNOW        (49UL)  Just abort everything
    % MDH_RETRO_LASTHEARTBEAT       (50UL)  This adc is from the last heartbeat (a dummy)
    % MDH_RETRO_DUMMYSCAN           (51UL)  This adc is just a dummy scan, throw it away
    % MDH_RETRO_ARRDETDISABLED      (52UL)  Disable all arrhythmia detection when this bit is found
    % MDH_B1_CONTROLLOOP            (53UL)  Marks the readout as to be used for B1 Control Loop
    % MDH_SKIP_ONLINE_PHASCOR       (54UL)  Marks scans not to be online phase corrected, even if online phase correction is switched on
    % MDH_SKIP_REGRIDDING           (55UL)  Marks scans not to be regridded, even if regridding is switched on
   
        EvalInfo = {'ACQEND', ...                    %  1 (0UL)
                    'RTFEEDBACK', ...                %  2 (1UL)
                    'HPFEEDBACK', ...                %  3 (2UL)
                    'ONLINE', ...                    %  4 (3UL)
                    'OFFLINE', ...                   %  5 (4UL)
                    'SYNCDATA', ...                  %  6 (5UL)
                    'Seven', ...                     %  7 (6UL)
                    'Eight', ...                     %  8 (7UL)
                    'LASTSCANINCONCAT', ...          %  9 (8UL)
                    'Ten', ...                       % 10 (9UL)
                    'RAWDATACORRECTION', ...         % 11 (10UL)
                    'LASTSCANINMEAS', ...            % 12 (11UL)
                    'SCANSCALEFACTOR', ...           % 13 (12UL)
                    '2NDHADAMARPULSE', ...           % 14 (13UL)
                    'REFPHASESTABSCAN', ...          % 15 (14UL)
                    'PHASESTABSCAN', ...             % 16 (15UL)
                    'D3FFT', ...                     % 17 (16UL)
                    'SIGNREV', ...                   % 18 (17UL)
                    'PHASEFFT', ...                  % 19 (18UL)
                    'SWAPPED', ...                   % 20 (19UL)
                    'POSTSHAREDLINE', ...            % 21 (20UL)
                    'PHASCOR', ...                   % 22 (21UL)
                    'PATREFSCAN', ...                % 23 (22UL)
                    'PATREFANDIMASCAN', ...          % 24 (23UL)
                    'REFLECT', ...                   % 25 (24UL)
                    'NOISEADJSCAN', ...              % 26 (25UL)
                    'SHARENOW', ...                  % 27 (26UL)
                    'LASTMEASUREDLINE', ...          % 28 (27UL)
                    'FIRSTSCANINSLICE', ...          % 29 (28UL)
                    'LASTSCANINSLICE', ...           % 30 (29UL)
                    'TREFFECTIVEBEGIN', ...          % 31 (30UL)
                    'TREFFECTIVEEND'};               % 32 (31UL)

        EvalInfo2 ={'MDS_REF_POSITION', ...          %  1 (32UL)
                    'SLC_AVERAGED', ...              %  2 (33UL)
                    'TAGFLAG1', ...                  %  3 (34UL)
                    'CT_NORMALIZE', ...              %  4 (35UL)
                    'SCAN_FIRST', ...                %  5 (36UL)
                    'SCAN_LAST', ...                 %  6 (37UL)
                    'Seven', ...                     %  7
                    'Eight', ...                     %  8
                    'FIRST_SCAN_IN_BLADE', ...       %  9 (40UL)
                    'LAST_SCAN_IN_BLADE', ...        % 10 (41UL)
                    'LAST_BLADE_IN_TR', ...          % 11 (42UL)
                    'Twelve', ...                    % 12
                    'PACE', ...                      % 13 (44UL)
                    'RETRO_LASTPHASE', ...           % 14 (45UL)
                    'RETRO_ENDOFMEAS', ...           % 15 (46UL)
                    'RETRO_REPEATTHISHEARTBEAT', ... % 16 (47UL)
                    'RETRO_REPEATPREVHEARTBEAT', ... % 17 (48UL)
                    'RETRO_ABORTSCANNOW', ...        % 18 (49UL)
                    'RETRO_LASTHEARTBEAT', ...       % 19 (50UL)
                    'RETRO_DUMMYSCAN', ...           % 20 (51UL)
                    'RETRO_ARRDETDISABLED', ...      % 21 (52UL)
                    'B1_CONTROLLOOP', ...            % 22 (53UL)
                    'SKIP_ONLINE_PHASCOR', ...       % 23 (54UL)
                    'SKIP_REGRIDDING'};              % 24 (55UL)
   
        MDH = '';
        for kb = 1:32
            if bitget(h.Mdh.EvalInfoMask, kb)
                MDH = [ MDH ' ' EvalInfo{kb} ];
            end
        end
        for kb = 1:24
            if bitget(h.Mdh.EvalInfoMask2, kb)
                MDH = [ MDH ' ' EvalInfo2{kb} ];
            end
        end
        
    end % of infoMDH

    function setWaitbar(h)
    %% select progress bar
        if h.WithWaitBar < 0, h.WithWaitBar = 0; end
    end
    
    function updateWaitbar(h)
    %% update progress bar
        if h.WithWaitBar < 0, return; end
        if h.BufferedDMA
            x = round( 100 * (h.File.pos + (h.BufferPos - h.BufferSiz)*4)  / h.File.size) / 100;
        else
            x = round( 100 * h.File.pos / h.File.size) / 100;
        end

        if isempty(h.WaitBar)
            h.WithWaitBar = x;
            if h.File.syngo == h.syngoVA
                [~,str] = fileparts(h.File.path);
                str = [ str ' (' h.File.SeriesDescription ')' ];
            else
                str = h.File.basename;
            end
            if length(str)>128, str = [ str(1:124) ' ...' ]; end
            h.WaitBar = waitbar(x, str, ...
                'CreateCancelBtn', {@stopit, h} ) ;
            set(findall(h.WaitBar,'type','text'),'Interpreter','none');
            pause(0.01);
        elseif x ~= h.WithWaitBar
            h.WithWaitBar = x;
            waitbar(x, h.WaitBar);
        end
    end % end of updateWaitbar
    
    function Progress = filePos(h)
        Progress = 100 * h.File.pos/ h.File.size;
        % if nargout == 0
        %    fprintf('%.2f %%, filepos = %d of %d bytes, buffer size = %d ulong\n', ...
        %        Progress, h.File.pos, h.File.size, h.BufferSiz);
        % end
    end
    
    function version(h)
        disp(h.Version);
    end
    
  end % of public methods in getwix class

  % /////////////////////////////////////////////////////////////////////////
  methods(Access = private)

    function Status = fetchBuffer(h, MnSizInLong)
    %% used by fetchReadout to store twix file into a buffer
        Status = true;

        if (h.BufferSiz - h.BufferPos + 1) < MnSizInLong
            
            if h.File.pos + h.BufferMxByte > h.File.size
                BufSiz = (h.File.size - h.File.pos) / 4;
                if BufSiz <= 0 % end of file
                    close(h);
                    Status = false;
                    return
                end
            else
                BufSiz = h.BufferMxByte / 4;
            end

            RestBuf = [];
            if h.BufferSiz > 0 && h.BufferPos <= h.BufferSiz
                RestBuf = h.Buffer(h.BufferPos:h.BufferSiz);
            end
            Buf   = fread(h.fid, BufSiz, 'ulong=>uint32');
            
            if feof(h.fid)
                Status = false;
                filePos(h);
                disp( [ mfilename('class') ': error while reading: ' fopen(h.fid) ] );
                return
            end
            
            h.File.pos  = h.File.pos + BufSiz * 4;

            if ~isempty(RestBuf), Buf = [ RestBuf; Buf ];  end
            h.Buffer    = Buf;
            h.BufferPos = 1;
            h.BufferSiz = length(h.Buffer);
            
            if h.WithWaitBar < 0, filePos(h); end
        end
        
    end
    
    function nCoils = getReadoutByRead(h)
    %% classical approach using direct file read access
        
        if h.File.syngo == h.syngoVD
            MdhSize = 192;
        else
            MdhSize = 128;
        end
        
        % convert each channels
        nCoils         = 0;
        nExpectedCoils = 1;
        kc             = 1;
        
        while (kc <= nExpectedCoils)
            
            % in VE11 only one measurement data header for all channels
            if kc == 1 || h.File.syngo ~= h.syngoVD
                
                % First, read the mdh (measurement data header)
                DMAlen1 = fread(h.fid, 1, 'ushort'); % this probably includes all channels
                if feof(h.fid), return; end
                DMAlen2 = fread(h.fid, 2, 'uint8');
                % fprintf(1, 'kc=%d DMA = %d  %d  %d ', kc, DMAlen1, DMAlen2(1), DMAlen2(2));
                
                h.Mdh.DMAlength               = DMAlen1 + DMAlen2(1) * 65536; % 256*256 = 65536
                h.Mdh.DMAflags                = DMAlen2(2);
                h.Mdh.MeasUID                 = fread(h.fid, 1, 'long');
                h.Mdh.ScanCounter             = fread(h.fid, 1, 'ulong');
                
                % time since 00:00 in 2.5 ms ticks
                h.Mdh.TimeStamp               = fread(h.fid, 1, 'ulong');
                
                % time since last trigger in 2.5 ms ticks
                h.Mdh.PMUTimeStamp            = fread(h.fid, 1, 'ulong');
                
                if h.File.syngo == h.syngoVD
                    % VE11 structure
                    h.Mdh.SystemType          = fread(h.fid, 1, 'ushort');
                    h.Mdh.PTABPos.Delay       = fread(h.fid, 1, 'ushort');
                    h.Mdh.PTABPos.X           = fread(h.fid, 1, 'long');
                    h.Mdh.PTABPos.Y           = fread(h.fid, 1, 'long');
                    h.Mdh.PTABPos.Z           = fread(h.fid, 1, 'long');
                    h.Mdh.Reserved	          = fread(h.fid, 1, 'ulong');
                end
                
                % EvalInfoMask
                h.Mdh.EvalInfoMask            = fread(h.fid, 1, 'ulong');
                
                % Last 32 bits of 64 bit EvalInfoMask
                h.Mdh.EvalInfoMask2           = fread(h.fid, 1, 'ulong');
                
                h.Mdh.SamplesInScan           = fread(h.fid, 1, 'ushort');
                h.Mdh.UsedChannels            = fread(h.fid, 1, 'ushort');
                
                h.Mdh.Loop.Clin               = fread(h.fid, 1, 'ushort'); % line
                h.Mdh.Loop.Cacq               = fread(h.fid, 1, 'ushort'); % average
                h.Mdh.Loop.Cslc               = fread(h.fid, 1, 'ushort'); % slice
                h.Mdh.Loop.Cpar               = fread(h.fid, 1, 'ushort'); % partition
                h.Mdh.Loop.Ceco               = fread(h.fid, 1, 'ushort'); % echo
                h.Mdh.Loop.Cphs               = fread(h.fid, 1, 'ushort'); % phase
                h.Mdh.Loop.Crep               = fread(h.fid, 1, 'ushort'); % repetition
                h.Mdh.Loop.Cset               = fread(h.fid, 1, 'ushort'); % set (TI/TD/Diff)
                h.Mdh.Loop.Cseg               = fread(h.fid, 1, 'ushort'); % segment (forward/reflected)
                h.Mdh.Loop.Cida               = fread(h.fid, 1, 'ushort'); % dimension A
                h.Mdh.Loop.Cidb               = fread(h.fid, 1, 'ushort'); % dimension B
                h.Mdh.Loop.Cidc               = fread(h.fid, 1, 'ushort'); % dimension C
                h.Mdh.Loop.Cidd               = fread(h.fid, 1, 'ushort'); % dimension D
                h.Mdh.Loop.Cide               = fread(h.fid, 1, 'ushort'); % dimension E
                
                h.Mdh.CutOffDataPre           = fread(h.fid, 1, 'ushort');
                h.Mdh.CutOffDataPost          = fread(h.fid, 1, 'ushort');
                
                h.Mdh.KSpaceCentreColumn      = fread(h.fid, 1, 'ushort');
                h.Mdh.CoilSelect              = fread(h.fid, 1, 'ushort');
                h.Mdh.ReadOutOffcentre        = fread(h.fid, 1, 'float');
                h.Mdh.TimeSinceLastRF         = fread(h.fid, 1, 'ulong');
                h.Mdh.KSpaceCentreLineNo      = fread(h.fid, 1, 'ushort');
                h.Mdh.KSpaceCentrePartitionNo = fread(h.fid, 1, 'ushort');
                
                if h.File.syngo == h.syngoVD
                    
                    % VE11 version - Slice data swapped with longer ICE params
                    h.Mdh.SD.SlicePos.Sag     = fread(h.fid, 1, 'float');
                    h.Mdh.SD.SlicePos.Cor     = fread(h.fid, 1, 'float');
                    h.Mdh.SD.SlicePos.Tra     = fread(h.fid, 1, 'float');
                    h.Mdh.SD.Quaternion       = fread(h.fid, 4, 'float');
                    h.Mdh.IceProgramPara      = fread(h.fid, 24, 'ushort');
                    h.Mdh.FreePara            = fread(h.fid, 4, 'ushort'); % reserved param
                    h.Mdh.ApplicationCounter  = fread(h.fid, 1, 'ushort');
                    h.Mdh.ApplicationMask     = fread(h.fid, 1, 'ushort');
                    h.Mdh.CRC			      = fread(h.fid, 1, 'ulong'); % check 32 sum
                    
                    % current position VE11: 192 bytes (MDH)
                    h.File.pos = h.File.pos + MdhSize;
                    
                else
                    
                    % VB15 version - ICE, slice data, channel ID, PTAB
                    h.Mdh.IceProgramPara      = fread(h.fid, 4, 'ushort');
                    h.Mdh.FreePara            = fread(h.fid, 4, 'ushort');
                    
                    h.Mdh.SD.SlicePos.Sag     = fread(h.fid, 1, 'float');
                    h.Mdh.SD.SlicePos.Cor     = fread(h.fid, 1, 'float');
                    h.Mdh.SD.SlicePos.Tra     = fread(h.fid, 1, 'float');
                    h.Mdh.SD.Quaternion       = fread(h.fid, 4, 'float');
                    
                    h.Cdh(kc).Id              = fread(h.fid, 1, 'ushort');
                    h.Mdh.PTABPosNeg          = fread(h.fid, 1, 'ushort');
                    
                    % current position VB17: 128 bytes (MDH)
                    h.File.pos = h.File.pos + MdhSize;
                    
                end
                
                if kc == 1
                    
                    % this allows to adapt the number of used channels
                    nExpectedCoils = h.Mdh.UsedChannels;
                    
                    % timestamp in ms to Mdh
                    h.Mdh.TimeStampInMs      = 2.5 * h.Mdh.TimeStamp;    % now in milliseconds
                    h.Mdh.PMUTimeStampInMs   = 2.5 * h.Mdh.PMUTimeStamp; % now in milliseconds
                    
                end
                
            end

            
            % skip MDH_SYNCDATA
            if MDH_SYNCDATA(h)
                nCoils = -1;
                fseek(h.fid, h.Mdh.DMAlength - MdhSize, 'cof');
                return
            end
                
            if h.File.syngo == h.syngoVD
                
                %% VE 11 - channel data header (Cdh)
                h.Cdh(kc).TypeLen          = fread(h.fid, 1, 'ulong');
                h.Cdh(kc).MeasUID          = fread(h.fid, 1, 'long');
                h.Cdh(kc).ScanCounter      = fread(h.fid, 1, 'ulong');
                h.Cdh(kc).Reserved         = fread(h.fid, 1, 'ulong');
                h.Cdh(kc).SequenceTime     = fread(h.fid, 1, 'ulong');
                h.Cdh(kc).Unused1          = fread(h.fid, 1, 'ulong');
                h.Cdh(kc).Id               = fread(h.fid, 1, 'ushort');
                h.Cdh(kc).Unused3          = fread(h.fid, 1, 'ushort');
                h.Cdh(kc).CRC              = fread(h.fid, 1, 'ulong');
                
                % current position VE11: 32 bytes (CDH)
                h.File.pos = h.File.pos + 32;
                
            end
            
            % collect ADC as complex float
            adc = fread(h.fid, 2 * h.Mdh.SamplesInScan, 'float');
            adc = reshape(adc, 2, h.Mdh.SamplesInScan);
            h.Cdh(kc).adc = complex(adc(1,:), adc(2,:));
            
            % current position VE11: ADC size (complex)
            h.File.pos = h.File.pos + 2 * h.Mdh.SamplesInScan * 4;
            
            nCoils = nCoils + 1;
            kc = kc + 1;
            
        end % loop channels
        
        
    end % end of classical direct file read access
    
    function nCoils = getReadoutByBuffer(h)
    %% buffered approach using type cast //////////////////////////////
        
        if h.File.syngo == h.syngoVD
            MdhSize = 192/4;
        else
            MdhSize = 128/4;
        end
        
        nCoils         = 0;
        nExpectedCoils = 1;
        kc             = 1;

        while (kc <= nExpectedCoils)

            if kc == 1 || h.File.syngo ~= h.syngoVD
                
                % check if buffer is empty and contains DMA length
                if ~fetchBuffer(h, MdhSize), return; end
                buf         = h.Buffer(h.BufferPos:(h.BufferPos+MdhSize-1));
                bush        = typecast( buf, 'uint16');
                h.BufferPos = h.BufferPos + MdhSize;
                
                if h.File.syngo == h.syngoVD
                    
                    % VE11 header is made of 192 bytes (48 unsigned longs)
                    
                    DMAlen1                       = double(bush(1));
                    DMAlen2                       = double(typecast( bush(2), 'uint8'));
                    h.Mdh.DMAlength               = DMAlen1 + DMAlen2(1) * 65536;
                    h.Mdh.DMAflags                = DMAlen2(2);

                    h.Mdh.MeasUID                 = typecast(buf(2), 'int32');
                    h.Mdh.ScanCounter             = double(buf(3));
                    h.Mdh.TimeStamp               = double(buf(4));
                    h.Mdh.PMUTimeStamp            = double(buf(5));
                    
                    h.Mdh.SystemType              = bush(11);
                    h.Mdh.PTABPos.Delay           = bush(12);
                    h.Mdh.PTABPos.X               = typecast(buf(7), 'int32');
                    h.Mdh.PTABPos.Y               = typecast(buf(8), 'int32');
                    h.Mdh.PTABPos.Z               = typecast(buf(9), 'int32');
                    h.Mdh.Reserved                = buf(10);
                      
                    h.Mdh.EvalInfoMask            = buf(11);
                    h.Mdh.EvalInfoMask2           = buf(12);
                    h.Mdh.SamplesInScan           = double(bush(25));
                    h.Mdh.UsedChannels            = double(bush(26));
                    h.Mdh.Loop.Clin               = double(bush(27));
                    h.Mdh.Loop.Cacq               = double(bush(28));
                    h.Mdh.Loop.Cslc               = double(bush(29));
                    h.Mdh.Loop.Cpar               = double(bush(30));
                    h.Mdh.Loop.Ceco               = double(bush(31));
                    h.Mdh.Loop.Cphs               = double(bush(32));
                    h.Mdh.Loop.Crep               = double(bush(33));
                    h.Mdh.Loop.Cset               = double(bush(34));
                    h.Mdh.Loop.Cseg               = double(bush(35));
                    h.Mdh.Loop.Cida               = double(bush(36));
                    h.Mdh.Loop.Cidb               = double(bush(37));
                    h.Mdh.Loop.Cidc               = double(bush(38));
                    h.Mdh.Loop.Cidd               = double(bush(39));
                    h.Mdh.Loop.Cide               = double(bush(40));
                    
                    h.Mdh.CutOffDataPre           = double(bush(41));
                    h.Mdh.CutOffDataPost          = double(bush(42));
                    h.Mdh.KSpaceCentreColumn      = double(bush(43));
                    h.Mdh.CoilSelect              = double(bush(44));
                    h.Mdh.ReadOutOffcentre        = double(typecast(buf(23), 'single'));
                    h.Mdh.TimeSinceLastRF         = buf(24);
                    h.Mdh.KSpaceCentreLineNo      = double(bush(49));
                    h.Mdh.KSpaceCentrePartitionNo = double(bush(50));
                    
                    h.Mdh.SD.SlicePos.Sag         = double(typecast(buf(26), 'single'));
                    h.Mdh.SD.SlicePos.Cor         = double(typecast(buf(27), 'single'));
                    h.Mdh.SD.SlicePos.Tra         = double(typecast(buf(28), 'single'));
                    h.Mdh.SD.Quaternion           = double(typecast(buf(29:32), 'single'));
                    
                    h.Mdh.IceProgramPara          = double(bush(65:88));
                    h.Mdh.FreePara                = double(bush(89:92));
                    
                    h.Mdh.ApplicationCounter      = double(bush(93));
                    h.Mdh.ApplicationMask         = double(bush(94));
                    h.Mdh.CRC                     = buf(48);

                else
                    
                    % VB17 header is made of 128 bytes (32 unsigned longs)
                    
                    DMAlen1                       = double(bush(1));
                    DMAlen2                       = double(typecast( bush(2), 'uint8'));
                    h.Mdh.DMAlength               = DMAlen1 + DMAlen2(1) * 65536;
                    h.Mdh.DMAflags                = DMAlen2(2);
                    
                    if kc == 1

                        h.Mdh.MeasUID             = typecast(buf(2), 'int32');
                        h.Mdh.ScanCounter         = double(buf(3));
                        h.Mdh.TimeStamp           = double(buf(4));
                        h.Mdh.PMUTimeStamp        = double(buf(5));
                        
                        h.Mdh.EvalInfoMask        = buf(6);
                        h.Mdh.EvalInfoMask2       = buf(7);
                        h.Mdh.SamplesInScan       = double(bush(15));
                        h.Mdh.UsedChannels        = double(bush(16));
                        h.Mdh.Loop.Clin           = double(bush(17));
                        h.Mdh.Loop.Cacq           = double(bush(18));
                        h.Mdh.Loop.Cslc           = double(bush(19));
                        h.Mdh.Loop.Cpar           = double(bush(20));
                        h.Mdh.Loop.Ceco           = double(bush(21));
                        h.Mdh.Loop.Cphs           = double(bush(22));
                        h.Mdh.Loop.Crep           = double(bush(23));
                        h.Mdh.Loop.Cset           = double(bush(24));
                        h.Mdh.Loop.Cseg           = double(bush(25));
                        h.Mdh.Loop.Cida           = double(bush(26));
                        h.Mdh.Loop.Cidb           = double(bush(27));
                        h.Mdh.Loop.Cidc           = double(bush(28));
                        h.Mdh.Loop.Cidd           = double(bush(29));
                        h.Mdh.Loop.Cide           = double(bush(30));
                        
                        h.Mdh.CutOffDataPre       = double(bush(31));
                        h.Mdh.CutOffDataPost      = double(bush(32));
                        h.Mdh.KSpaceCentreColumn  = double(bush(33));
                        h.Mdh.CoilSelect          = double(bush(34));
                        h.Mdh.ReadOutOffcentre    = double(typecast(buf(18), 'single'));
                        h.Mdh.TimeSinceLastRF     = double(buf(19));
                        h.Mdh.KSpaceCentreLineNo  = double(bush(39));
                        h.Mdh.KSpaceCentrePartitionNo = double(bush(40));
                        
                        h.Mdh.IceProgramPara      = double(bush(41:44));
                        h.Mdh.FreePara            = double(bush(45:48));
                        
                        h.Mdh.SD.SlicePos.Sag     = double(typecast(buf(25), 'single'));
                        h.Mdh.SD.SlicePos.Cor     = double(typecast(buf(26), 'single'));
                        h.Mdh.SD.SlicePos.Tra     = double(typecast(buf(27), 'single'));
                        h.Mdh.SD.Quaternion       = double(typecast(buf(28:31),'single'));
                        h.Mdh.PTABPosNeg          = double(bush(64));

                    end
                    % ChannelId is copied to channel header
                    h.Cdh(kc).ChannelId           = bush(63);

                end
                
                % this allows to adapt the number of used channels
                nExpectedCoils = h.Mdh.UsedChannels;                
                
                % timestamp in ms to Mdh
                h.Mdh.TimeStampInMs      = 2.5 * h.Mdh.TimeStamp;    % now in milliseconds
                h.Mdh.PMUTimeStampInMs   = 2.5 * h.Mdh.PMUTimeStamp; % now in milliseconds
                
            end
            
            % skip MDH_SYNCDATA
            if MDH_SYNCDATA(h), return; end
            
            % ADC is made of 4 bytes * h.Mdh.SamplesInScan * 2
            nAdc        = 2 * h.Mdh.SamplesInScan;
            nBuf        = nAdc;
            if h.File.syngo == h.syngoVD, nBuf = nBuf + 8; end
           
            if ~fetchBuffer(h, nBuf), return; end
            ipos        = h.BufferPos;
            h.BufferPos = h.BufferPos + nBuf; % update buffer position
        
            if h.File.syngo == h.syngoVD

                buf     = h.Buffer(ipos:(ipos+7));
                bush    = typecast( buf, 'uint16');
                               
                %% VE 11 - channel data header (Cdh)
                h.Cdh(kc).TypeLen        = buf(1);
                h.Cdh(kc).MeasUID        = typecast(buf(2), 'int32');
                h.Cdh(kc).ScanCounter    = buf(3);
                h.Cdh(kc).Reserved       = buf(4);
                h.Cdh(kc).SequenceTime   = buf(5);
                h.Cdh(kc).Unused1        = buf(6);
                h.Cdh(kc).Id             = bush(14);
                h.Cdh(kc).Unused3        = bush(15);
                h.Cdh(kc).CRC            = buf(8);
                
                ipos    = ipos + 8;
                
            end
            
            % collect ADC as complex float
            adc   = double( reshape( ...
                      typecast(h.Buffer(ipos:(ipos+nAdc-1)), 'single'), ...
                      2, h.Mdh.SamplesInScan) );
            h.Cdh(kc).adc = complex(adc(1,:), adc(2,:));
            
            nCoils      = nCoils + 1;
            kc          = kc + 1;
        end
        
    end
   
    
    function Status = fetchMrProt(h, hdr)
    %% fetchMrProt(h, hdr) -- extract MrProtocol from TWIX

        Status = false;
        h.File.FieldNames = {};
        nFields = 0;

        ibeg = strfind(hdr, ['### ASCCONV BEGIN ###' char(10) 'ulVersion'] );
        % guess VE11
        if isempty(ibeg)
            ibeg = strfind(hdr, '### ASCCONV BEGIN object=MrProt' );
        end
        % guess VA25
        if isempty(ibeg)
            ibeg = strfind(hdr, '### ASCCONV BEGIN ###' );
        end
        if isempty(ibeg)
            h.File.Msg = [ mfilename('class') ': MrProtocol not found' ];
            return
        end

        iend = strfind(hdr(ibeg(1):end),'### ASCCONV END ###');
        if isempty(iend)
            iend = length(hdr);
        else
            iend = ibeg(1) + iend + 17;
        end

        pro = [];
        
        hdr = hdr(ibeg(1):iend(1));
        hdr(length(hdr)+1) = char(10);
        idx = strfind(hdr,char(10));
        istart = 1;
        for k=1:length(idx)

            str = hdr(istart:idx(k)-1);
            istart = idx(k)+1;
            ieq = strfind(str, ' = ');
            if length(ieq) == 1
                key  = deblank(str(1:ieq));
                val  = deblank(str(ieq+3:end));

                if isempty( strfind(key, '__attribute__') ) % guess VE11

                    if length(val)>1 && val(1) == char(9), val = val(2:end); end
                    key = strrep(key, '__', '');

                    iBra = strfind(key,'[');
                    iKet = strfind(key,']');
                    while (~isempty(iBra) && ~isempty(iKet))
                        key = [ key(1:iBra(1)-1) '(' ...
                            num2str( str2num( key(iBra(1)+1:iKet(1)-1) ) + 1 ) ...
                            ')' key(iKet(1)+1:end) ]; %#ok<ST2NM>
                        iBra = strfind(key,'[');
                        iKet = strfind(key,']');
                    end
                    if val(1)=='"' && val(end)=='"'
                        val = [ '''' val(2:end-1) '''' ];
                    end
                    if length(val)>2 && val(1)=='0' && val(2)=='x'
                        try
                            val = num2str(hex2dec(val(3:end)));
                        catch
                            disp([ mfilename('class') ': ' key ' = ' val ]);
                        end
                    end
                    try 
                        eval(['pro.' key ' = ' val ';']);
                        nFields = nFields + 1;
                        h.File.FieldNames{nFields,1} = key;
                    catch
                        disp([ mfilename('class') ': failed: ' key ' = ' val ]);
                    end
                end
            end
        end

        % add some zero default values, complete slice objects
        if ~isfield(pro,'lRepetitions'), pro.lRepetitions = 0; end
        for ks = 1:pro.sSliceArray.lSize
            asSlice(ks) = sliceData(h, pro.sSliceArray.asSlice(ks) );
        end
        pro.sSliceArray.asSlice = asSlice;
        for ks = 1:pro.sGroupArray.lSize
            asGroup(ks) = pro.sGroupArray.asGroup(ks);
            if ~isfield(asGroup(ks),'dDistFact')
                asGroup(ks).dDistFact = 0;
            end
        end
        pro.sGroupArray.asGroup = asGroup;
        
        h.MrProt = pro;
        
        % populate DICOM some fields
        h.File.SeriesDescription = MrProtString(h, pro.tProtocolName); % (0008,103e) LO

        try
            switch h.MrProt.sSliceArray.ucMode
                case 1, h.File.sliceOrderInfo    = 'ascending';
                case 2, h.File.sliceOrderInfo    = 'descending';
                case 4, h.File.sliceOrderInfo    = 'interleaved';
                case 8, h.File.sliceOrderInfo    = 'automatic';
                case 16, h.File.sliceOrderInfo   = 'apextobase';
                case 32, h.File.sliceOrderInfo   = 'basetoapex';
                otherwise, h.File.sliceOrderInfo = 'unknown';
            end
        catch
            h.File.sliceOrderInfo = 'unknown';
        end
 
        Status = true;
    end % of fetchMrProt

    function SL = sliceData(~, SL)
    %% complete dTra, dCor and dSag in slice data
        if ~isfield(SL,'sPosition'), SL.sPosition = []; end
        if ~isfield(SL.sPosition,'dTra'), SL.sPosition.dTra = 0; end
        if ~isfield(SL.sPosition,'dCor'), SL.sPosition.dCor = 0; end
        if ~isfield(SL.sPosition,'dSag'), SL.sPosition.dSag = 0; end
        if ~isfield(SL,'sNormal'), SL.sNormal = []; end
        if ~isfield(SL.sNormal,'dTra'), SL.sNormal.dTra = 0; end
        if ~isfield(SL.sNormal,'dCor'), SL.sNormal.dCor = 0; end
        if ~isfield(SL.sNormal,'dSag'), SL.sNormal.dSag = 0; end
        if ~isfield(SL,'dThickness'), SL.dThickness = 0; end
        if ~isfield(SL,'dPhaseFOV'), SL.dPhaseFOV = 0; end
        if ~isfield(SL,'dReadoutFOV'), SL.dReadoutFOV = 0; end
        if ~isfield(SL,'dInPlaneRot'), SL.dInPlaneRot = 0; end
    end

    function st = MrProtString(~,st)
    %% MrProtString(h, str) -- parse Siemens MrProtocol string for display
        if strcmp(st(1:2),'""'), st=st(3:end); end
        if strcmp(st(1),'"'), st=st(2:end); end
        if strcmp(st(end-1:end),'""'), st=st(1:end-2); end
        if strcmp(st(end:end),'"'), st=st(1:end-1); end
        v = regexp(st, '+AF8-', 'split');
        st=[v{1}]; for i=2:length(v), st=[st '_' v{i}]; end %#ok<*AGROW>
        v = regexp(st, '+AEA-', 'split');
        st=[v{1}]; for i=2:length(v), st=[st '@' v{i}]; end
    end % of MrProtString

    function fetchMrXML(h, hdr)
    %% fetchMrXML(h, header) -- parse XProtocol in XML

        % fp = fopen('XProtParsed.txt','w');
        nFields = length( h.File.FieldNames );
    
        x = [];
        
        idx1= strfind(hdr,'<Pipe');
        idx2= strfind(hdr,'<Param');
        idx=sort([idx1 idx2]);
        step.level = 0;
        step.val   = 'x';
        curidx     = 1;
        
        curlevel   = length( strfind( hdr(1:idx(1)+1), '{') ) - ...
                     length( strfind( hdr(1:idx(1)+1), '}') );
    
        for k=1:length(idx)-1
            
            k0  = idx(k)+1;
            seg = hdr(k0 : idx(k+1)-1);
            kk  = strfind(seg,'."');
            kkk = strfind(seg,'">');
            
            if ~isempty(kk) && ~isempty(kkk) && kk(1)+1 < kkk(1)
                P = seg(1:kk(1)-1);
                N = seg(kk(1)+2:kkk(1)-1);
                
                % compute difference in level
                level    = curlevel;
                curlevel = curlevel + ...
                           length( strfind(seg,'{') ) - ...
                           length( strfind(seg,'}') );
                    
                if  ~isempty(N)

                    % fprintf(fp,' P=%s N=%s level=%d curidx=%d -> level=%d val=%s k0=%d\n', ...
                    %            P,N,level,curidx,step(curidx).level,step(curidx).val,k0);
                    
                    for kk=curidx:-1:2
                        if level <= step(kk).level, curidx = kk-1; end
                    end

                    if  ( strcmp(P, 'ParamMap') || strcmp(P, 'ParamFunctor') || strcmp(P, 'ParamArray') || ...
                          strcmp(P, 'PipeService') || strcmp(P, 'Pipe') );

                        NN = [ step(curidx).val '.' N ];
                        curidx = curidx + 1;
                        step(curidx).val = NN;
                        step(curidx).level = level;

                    else

                        st = '';

                        % respect field name convention
                        if N(1) == '_' || isstrprop(N(1), 'digit')
                            N = ['X' N];
                        end

                        if strcmp(P, 'ParamString')
                            [~,val]=strtok(seg,'{');
                            val=strtok(val(2:end),'}');
                            [~,val]=strtok(val,'"'); %#ok<STTOK>
                            val=strtok(val,'"');
                            st = sprintf('%s.%s = ''%s'';', step(curidx).val, N, val);
                        elseif strcmp(P, 'ParamLong')
                            [~,val]=strtok(seg,'{');
                            val=strtok(val(2:end),'}');
                            val = strrep(val, sprintf('\n'), ' ');
                            st = sprintf('%s.%s = [%s];',step(curidx).val, N, val);
                        elseif strcmp(P, 'ParamDouble')
                            [~,val]=strtok(seg,'{');
                            val=strtok(val(2:end),'}');
                            val = strrep(val, sprintf('\n'), ' ');

                            kkk = strfind(val,'<Precision>');
                            if ~isempty(kkk)
                                val = val(kkk(1)+11:end);
                                [~, val]=strtok(val); %#ok<STTOK>
                            end
                            st = sprintf('%s.%s = [%s];',step(curidx).val, N, val);
                        elseif strcmp(P, 'ParamBool')
                            [~,val]=strtok(seg,'{');
                            val=strtok(val(2:end),'}');
                            [~,val]=strtok(val,'"'); %#ok<STTOK>
                            val=strtok(val,'"');
                            st = sprintf('%s.%s = [%s];',step(curidx).val, N, val);
                        end
                        if ~isempty(st)
                            try
                                eval(st)
                                nFields = nFields + 1;
                                h.File.FieldNames{nFields,1}=sprintf('XProt.%s.%s',step(curidx).val(3:end), N);
                                % fprintf(fp,'%s\n',st);
                            catch
                                % disp([ mfilename('class') ':XProt failed: ' st]);
                            end
                        end
                    end
                end
                
            end
        end
        h.XProt = x;
        % fclose(fp);

        try
            val = x.YAPS.tFrameOfReference;
            kdx = strfind(val,'.');
            val = val( kdx(10)+1 : kdx(11)-1 );
            h.File.SeriesDate = val(1:8);
            h.File.SeriesTime = [ val(9:14) '.' val(15:16) ];
        catch
            h.File.SeriesDate = '';
            h.File.SeriesTime = '';
        end
    end % of fetchMrXML
   
    function delete(h)
    %% destructor of getwix class
        if h.fid > -1
            disp([ mfilename('class') '::destructor closing ' fopen(h.fid) ])
            fclose(h.fid);
        end
    end % of getwix destructor

  end % end of private methods

end % end of classdef getwix

%% call back function
function stopit(~, ~, h)
    disp([ 'Interupting ' h.File.basename ])
    h.close(true)
end

function readTwix(h)
%% read() -- classic read of Siemens NUMARIS 4 twix files
% See also fetchReadout

if h.fid < 0
    h.Info = [];
    disp( [ mfilename('class') ': no current twix file' ] );
    return;
end

nCoils = h.fetchReadout();

h.pub.Info = struct( 'nCha', h.Mdh.UsedChannels, ...
    'nCol', h.Mdh.SamplesInScan, 'nLin', 0, 'nAcq', 0, 'nSlc', 0,...
    'nPar', 0, 'nEco', 0, 'nPhs', 0, 'nRep', 0, 'nSet', 0, 'nSeg', 0, ...
    'nIda', 0, 'nIdb', 0, 'nIdc', 0, 'nIdd', 0, 'nIde', 0);

nMdh = 0;
while nCoils > 0

    if h.pub.Info.nLin <= h.lin, h.pub.Info.nLin = h.lin; end
    if h.pub.Info.nAcq <= h.acq, h.pub.Info.nAcq = h.acq; end
    if h.pub.Info.nSlc <= h.slc, h.pub.Info.nSlc = h.slc; end
    if h.pub.Info.nPar <= h.par, h.pub.Info.nPar = h.par; end
    if h.pub.Info.nEco <= h.eco, h.pub.Info.nEco = h.eco; end
    if h.pub.Info.nPhs <= h.phs, h.pub.Info.nPhs = h.phs; end
    if h.pub.Info.nRep <= h.rep, h.pub.Info.nRep = h.rep; end
    if h.pub.Info.nSet <= h.set, h.pub.Info.nSet = h.set; end
    if h.pub.Info.nSeg <= h.seg, h.pub.Info.nSeg = h.seg; end
    if h.pub.Info.nIda <= h.ida, h.pub.Info.nIda = h.ida; end
    if h.pub.Info.nIdb <= h.idb, h.pub.Info.nIdb = h.idb; end
    if h.pub.Info.nIdc <= h.idc, h.pub.Info.nIdc = h.idc; end
    if h.pub.Info.nIdd <= h.idd, h.pub.Info.nIdd = h.idd; end
    if h.pub.Info.nIde <= h.ide, h.pub.Info.nIde = h.ide; end
    
    nMdh = nMdh + 1;
    h.pub.DMA(nMdh).Mdh = h.Mdh;
    h.pub.DMA(nMdh).Cdh = h.Cdh;
    nCoils = h.fetchReadout();
end

disp( [ mfilename('class') ':read in h.pub.DMA ' num2str(nMdh) ' Mdh+ADC from ' fopen(h.fid) ] );
h.close();
end % of read



%% help functions
function viewTwix(h)

    nCoils = h.fetchReadout();
   
    % skip SYNCDATA
    while nCoils > 0 && ~h.MDH_ACQEND && h.MDH_SYNCDATA
        nCoils = h.fetchReadout();
    end
    
    % guess k-space data size
    nScans = floor((h.File.size - h.File.oskip) / h.Mdh.DMAlength);

    h.pub.ksp = zeros(nScans, nCoils * h.Mdh.SamplesInScan);
    
    % open display
    fg = figure('name', [mfilename('class') '.view']);
    imagesc(log(h.pub.ksp));
    if h.File.syngo == h.syngoVA
        [~,str] = fileparts(h.File.path);
        str = [ str '/' h.File.basename ' (' h.File.SeriesDescription ')' ];
    else
        str = h.File.basename;
    end
    title  (str, 'Interpreter', 'none')
    if nCoils == 1
        xlabel ('ADC')
    else
        xlabel (sprintf('ADC(1) ... ADC(%d)', nCoils))
    end
    ylabel ('ScanCounter')
    pause(0.05)
    hold on

    % loop
    while nCoils > 0 && ~h.MDH_ACQEND
        for kc = 1:nCoils
            h.pub.ksp(h.Mdh.ScanCounter, ...
               ((kc-1)*h.Mdh.SamplesInScan+1):(kc*h.Mdh.SamplesInScan)) ...
               = abs(h.Cdh(kc).adc);
        end
        nCoils = h.fetchReadout();
        if mod(h.Mdh.ScanCounter,512)==1
            imagesc(log(h.pub.ksp));
            pause(0.05)
        end
    end
    figure(fg)
    imagesc(log(h.pub.ksp));
    hold off
    h.close();
end


 
%% map twix file into one array
function mapTwix(h)

    nCoils = h.fetchReadout();
   
    % skip SYNCDATA
    while nCoils > 0 && ~h.MDH_ACQEND && h.MDH_SYNCDATA
        nCoils = h.fetchReadout();
    end

    % prepare map dimensions
    h.pub.Info.nCha = h.Mdh.UsedChannels;
    h.pub.Info.nCol = h.Mdh.SamplesInScan;
    h.pub.Info.nLin = h.lin;
    h.pub.Info.nAcq = h.acq;
    h.pub.Info.nSlc = h.slc;
    h.pub.Info.nPar = h.par;
    h.pub.Info.nEco = h.eco;
    h.pub.Info.nPhs = h.phs;
    h.pub.Info.nRep = h.rep;
    h.pub.Info.nSet = h.set;
    h.pub.Info.nSeg = h.seg;
    h.pub.Info.nIda = h.ida;
    h.pub.Info.nIdb = h.idb;
    h.pub.Info.nIdc = h.idc;
    h.pub.Info.nIdd = h.idd;
    h.pub.Info.nIde = h.ide;

    % read map dimensions from XProt
    try
        h.getXProt();
        FBR = h.XProt.PARC.PIPE.acquisition.FeedbackRoot;
        if isfield(FBR,'NLinMeas') && ~isempty(FBR.NLinMeas), h.pub.Info.nLin = FBR.NLinMeas; end
        if isfield(FBR,'NAveMeas') && ~isempty(FBR.NAveMeas), h.pub.Info.nAcq = FBR.NAveMeas; end
        if isfield(FBR,'NSlcMeas') && ~isempty(FBR.NSlcMeas), h.pub.Info.nSlc = FBR.NSlcMeas; end
        if isfield(FBR,'NParMeas') && ~isempty(FBR.NParMeas), h.pub.Info.nPar = FBR.NParMeas; end
        if isfield(FBR,'NEcoMeas') && ~isempty(FBR.NEcoMeas), h.pub.Info.nEco = FBR.NEcoMeas; end
        if isfield(FBR,'NPhsMeas') && ~isempty(FBR.NPhsMeas), h.pub.Info.nPhs = FBR.NPhsMeas; end
        if isfield(FBR,'NRepMeas') && ~isempty(FBR.NRepMeas), h.pub.Info.nRep = FBR.NRepMeas; end
        if isfield(FBR,'NSetMeas') && ~isempty(FBR.NSetMeas), h.pub.Info.nSet = FBR.NSetMeas; end
        if isfield(FBR,'NSegMeas') && ~isempty(FBR.NSegMeas), h.pub.Info.nSeg = FBR.NSegMeas; end
        if isfield(FBR,'NIdaMeas') && ~isempty(FBR.NIdaMeas), h.pub.Info.nIda = FBR.NIdaMeas; end
        if isfield(FBR,'NIdbMeas') && ~isempty(FBR.NIdbMeas), h.pub.Info.nIdb = FBR.NIdbMeas; end
        if isfield(FBR,'NIdcMeas') && ~isempty(FBR.NIdcMeas), h.pub.Info.nIdc = FBR.NIdcMeas; end
        if isfield(FBR,'NIddMeas') && ~isempty(FBR.NIddMeas), h.pub.Info.nIdd = FBR.NIddMeas; end
        if isfield(FBR,'NIdeMeas') && ~isempty(FBR.NIdeMeas), h.pub.Info.nIde = FBR.NIdeMeas; end
    catch
    end

    % declare map
    h.pub.map = zeros(h.pub.Info.nCha, ...
                      h.pub.Info.nCol, ...
                      h.pub.Info.nLin, ...
                      h.pub.Info.nAcq, ...
                      h.pub.Info.nSlc, ...
                      h.pub.Info.nPar, ...
                      h.pub.Info.nEco, ...
                      h.pub.Info.nPhs, ...
                      h.pub.Info.nRep, ...
                      h.pub.Info.nSet, ...
                      h.pub.Info.nSeg, ...
                      h.pub.Info.nIda, ...
                      h.pub.Info.nIdb, ...
                      h.pub.Info.nIdc, ...
                      h.pub.Info.nIdd, ...
                      h.pub.Info.nIde);
                  
    % declare evalInfoMask                  
    h.pub.evalInfoMask = zeros(h.pub.Info.nLin, ...
                      h.pub.Info.nAcq, ...
                      h.pub.Info.nSlc, ...
                      h.pub.Info.nPar, ...
                      h.pub.Info.nEco, ...
                      h.pub.Info.nPhs, ...
                      h.pub.Info.nRep, ...
                      h.pub.Info.nSet, ...
                      h.pub.Info.nSeg, ...
                      h.pub.Info.nIda, ...
                      h.pub.Info.nIdb, ...
                      h.pub.Info.nIdc, ...
                      h.pub.Info.nIdd, ...
                      h.pub.Info.nIde, ...
                      'uint32');
  
    % loop
    while nCoils > 0 && ~h.MDH_ACQEND
        h.pub.evalInfoMask(h.lin, ...
                           h.acq, ...
                           h.slc, ...
                           h.par, ...
                           h.eco, ...
                           h.phs, ...
                           h.rep, ...
                           h.set, ...
                           h.seg, ...
                           h.ida, ...
                           h.idb, ...
                           h.idc, ...
                           h.idd, ...
                           h.ide) = h.Mdh.EvalInfoMask;
        for kc = 1:nCoils
            h.pub.map(kc, ...
                      :, ...
                      h.lin, ...
                      h.acq, ...
                      h.slc, ...
                      h.par, ...
                      h.eco, ...
                      h.phs, ...
                      h.rep, ...
                      h.set, ...
                      h.seg, ...
                      h.ida, ...
                      h.idb, ...
                      h.idc, ...
                      h.idd, ...
                      h.ide) = h.Cdh(kc).adc(:);
        end
        nCoils = h.fetchReadout();
    end
    h.close();

    % update loop dimensions
    dim             = size(h.pub.map);  dim(17) = 0;
    h.pub.Info.nCha = dim(1);
    h.pub.Info.nCol = dim(2);
    h.pub.Info.nLin = dim(3);
    h.pub.Info.nAcq = dim(4);
    h.pub.Info.nSlc = dim(5);
    h.pub.Info.nPar = dim(6);
    h.pub.Info.nEco = dim(7);
    h.pub.Info.nPhs = dim(8);
    h.pub.Info.nRep = dim(9);
    h.pub.Info.nSet = dim(10);
    h.pub.Info.nSeg = dim(11);
    h.pub.Info.nIda = dim(12);
    h.pub.Info.nIdb = dim(13);
    h.pub.Info.nIdc = dim(14);
    h.pub.Info.nIdd = dim(15);
    h.pub.Info.nIde = dim(16);
    
    % prepare map reshape
    str = { 'Cha', 'Col' };
    dim = [ h.pub.Info.nCha, h.pub.Info.nCol ];
    n   = 2;
    
    % prune dimensions equal to zero or one
    if h.pub.Info.nLin>1, n=n+1; str{n}='Lin'; dim(n)=h.pub.Info.nLin; end
    if h.pub.Info.nAcq>1, n=n+1; str{n}='Acq'; dim(n)=h.pub.Info.nAcq; end
    if h.pub.Info.nSlc>1, n=n+1; str{n}='Slc'; dim(n)=h.pub.Info.nSlc; end
    if h.pub.Info.nPar>1, n=n+1; str{n}='Par'; dim(n)=h.pub.Info.nPar; end
    if h.pub.Info.nEco>1, n=n+1; str{n}='Eco'; dim(n)=h.pub.Info.nEco; end
    if h.pub.Info.nPhs>1, n=n+1; str{n}='Phs'; dim(n)=h.pub.Info.nPhs; end
    if h.pub.Info.nRep>1, n=n+1; str{n}='Rep'; dim(n)=h.pub.Info.nRep; end
    if h.pub.Info.nSet>1, n=n+1; str{n}='Set'; dim(n)=h.pub.Info.nSet; end
    if h.pub.Info.nSeg>1, n=n+1; str{n}='Seg'; dim(n)=h.pub.Info.nSeg; end
    if h.pub.Info.nIda>1, n=n+1; str{n}='Ida'; dim(n)=h.pub.Info.nIda; end
    if h.pub.Info.nIdb>1, n=n+1; str{n}='Idb'; dim(n)=h.pub.Info.nIdb; end
    if h.pub.Info.nIdc>1, n=n+1; str{n}='Idc'; dim(n)=h.pub.Info.nIdc; end
    if h.pub.Info.nIdd>1, n=n+1; str{n}='Idd'; dim(n)=h.pub.Info.nIdd; end
    if h.pub.Info.nIde>1, n=n+1; str{n}='Ide'; dim(n)=h.pub.Info.nIde; end

    % reshape twix map array
    h.pub.map          = reshape(h.pub.map, dim);
    
    % reshape EvalInfoMask array
    h.pub.evalInfoMask = reshape(h.pub.evalInfoMask, dim(3:end));
    
    h.pub.mapDim       = dim;
    h.pub.mapStr       = str;
    
    fprintf(1, '%s\nsize(h.pub.map):\n', h.File.basename);
    fprintf(1, '  %s',str{:});
    fprintf(1, '\n');
    fprintf(1,'%5d',dim);

    fprintf(1, '\nsize(h.pub.evalInfoMask):\n');
    fprintf(1, '  %s',str{3:end});
    fprintf(1, '\n');
    fprintf(1,'%5d',dim(3:end));
    fprintf(1, '\n');
end
