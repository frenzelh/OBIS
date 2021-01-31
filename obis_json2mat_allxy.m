% obis_json2mat_allxy - Processes OBIS files in JSON format:
% Binning to 1x1 degree WOA grid with 33 pre-WOD13 depth levels, 
% matching with T/S/O2 from WOA. Output is written to mat files.
% NOTE: Contains hard-coded paths.
%
% Synopsis:
%
% Inputs:
%
% Outputs:
%
% Author: H. Frenzel, School of Oceanography, UW
%
% Latest revision: October 15, 2019 (updates for OBIS API v3)
%
% First version: January 16, 2019

clear all;

% set defaults
in_dir = './JSON/';
out_dir = './OBIS_WOA_XY/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: READ WOA DATA (applied to all species)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WOA grid
load('/rraid1/frenzel/OBIS/woa_grid.mat','grid');
% only the grid vectors are needed
woa_grid.xt = grid.xt;
woa_grid.yt = grid.yt;
woa_grid.zt = grid.zw; % not grid.zt! those are lower box boundaries
clear grid; % uses over 500 MB!


% monthly WOA data
ncload('/rraid1/DATA/WOA18/woa18_monthly_33zlev_0-360.nc','temp','salt','o2')
temp(temp > 1e10) = nan;
salt(salt > 1e10) = nan;
o2(o2 > 1e10) = nan;
[NT,NZ,NY,NX] = size(o2);
WOA_mon.temp = temp; % tzyx
WOA_mon.salt = salt;
WOA_mon.o2 = o2;
clear temp salt o2;

% ETOPO topography (1 minute resolution)
fn_etopo = '/graid1/DATA/Topography/etopo1.cdf';
ncload(fn_etopo);

% delete the duplicated lon value at the ends and the corresponding topo
% values, as they fall into the middle of the circshifted fields
lon(end) = [];
topo(:,end) = [];

% put etopo fields into a structure for easy passing to function, adjust to
% lon = 0..360 degrees in the process
idx_neg = find(lon < 0);
lon(lon < 0) = lon(lon < 0) + 360.;
etopo.lon = circshift(lon, -length(idx_neg)); clear lon;

etopo.lat = lat; clear lat;
etopo.topo = circshift(topo, -length(idx_neg), 2); clear topo;

% now add back an end point so that every possible value of lon is within
% half the resolution of an etopo point
if (etopo.lon(1) > 1e-3)
    % should be 1/60
    etopo.lon = [etopo.lon(end) - 360; etopo.lon];
    etopo.topo = cat(2, etopo.topo(:,end), etopo.topo);
else
    warning('unexpected')
    keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2: READ OBIS DATA (loop over species)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files = dir([in_dir, '*_1.json']); 
nfiles = length(files);
first_file = 'Beroe_cucumis_1.json';
% numbers represented as strings
anc_fields_string = {'individualCount'};
% numbers represented as numbers
anc_fields_num = {'minimumDepthInMeters'; 'maximumDepthInMeters'};
% pure strings
anc_fields_char = {'occurrenceStatus'};
% all of these are eventually handled the same 
anc_fields_multi = [anc_fields_string; anc_fields_num; anc_fields_char];
% "species" may differ from file name!
anc_string = {'species';'family';'class';'phylum';'scientificNameID';...
    'scientificName';'originalScientificName';'scientificNameAuthorship'};
anc_num = {'aphiaID'; 'speciesID'}; % 'obisID'} % the last two should be the same
anc_fields_single = [anc_string; anc_num];
hw = waitbar(0,'');

if (exist('first_file', 'var') ~= 1 || isempty(first_file))
    process_files = 1;
else
    process_files = 0; % until first_file is encountered
end
for f = 1:nfiles
    fname = files(f).name
    if (~process_files)
        if (strcmp(fname, first_file))
            process_files = 1;
        else
            continue; % skip to next file
        end
    end
    msg = sprintf('Processing %s (%d/%d)', strrep(fname, '_', '\_'), ...
        f, nfiles);
    waitbar((f-1)/nfiles,hw,msg);
    index = 1;
    all_lon = [];
    all_lat = [];
    all_depth = [];
    all_year = [];
    all_month = [];
    all_temp = [];
    all_salt = [];
    all_o2 = [];
    all_z = [];
    all_y = [];
    all_x = [];
    all_zmin = [];
    all_zmax = [];
    all_bathy = [];
    is_first = 1; % used to assign "all_" char cells
    for i = 1:length(anc_fields_multi)
        eval(['all_',anc_fields_multi{i}, ' = [];']);
    end % i
    % ancillary data that are the same for all data points,
    % initialize only once (even if there are multiple files)
    for i = 1:length(anc_string)
        eval([anc_string{i}, ' = '''';']);
    end % i
    for i = 1:length(anc_num)
        eval([anc_num{i}, ' = nan;']);
    end % i
    format_counter = [0,0]; % D/M/Y or M/D/Y
    while (exist([in_dir, fname], 'file') == 2)
        json_data = jsondecode(fileread([in_dir, fname]));
        ndata = length(json_data);
        % core data
        lon = nan * ones(ndata,1);
        lat = nan * ones(ndata,1);
        depth = nan * ones(ndata,1);
        minimumDepthInMeters = nan * ones(ndata,1);
        maximumDepthInMeters = nan * ones(ndata,1);
        unix_year = nan * ones(ndata,1);
        year = nan * ones(ndata,1);
        month = nan * ones(ndata,1);
        % ancillary data that are different for data points
        eventDate  = strings(ndata,1); % treated differently below
        for i = 1:length(anc_fields_string)
            eval([anc_fields_string{i}, ' = nan * ones(ndata,1);']);
        end % i
        for i = 1:length(anc_fields_num)
            eval([anc_fields_num{i}, ' = nan * ones(ndata,1);']);
        end % i
        for i = 1:length(anc_fields_char)
            for d = 1:ndata
                eval([anc_fields_char{i}, '{d} = [];']);
            end % d
        end % i
        % loop over all data points from this file 
        for d = 1:ndata
            if (iscell(json_data))
                these_data = json_data{d};
            else
                these_data = json_data(d);
            end
            if (isfield(these_data, 'decimalLongitude'))
                lon(d) = these_data.decimalLongitude;
            end
            if (isfield(these_data, 'decimalLatitude'))
                lat(d) = these_data.decimalLatitude;
            end
            if (isfield(these_data, 'depth'))
                depth(d) = these_data.depth;
            end
            % "mid_date" is easily converted to y/m/d, but not all 
            % data points have the field
            if (isfield(these_data, 'date_mid'))
                unixdate=these_data.date_mid*1e-3; % convert millisec to sec
                dt = datetime(unixdate, 'ConvertFrom', 'posixtime');
                [year(d),month(d)] = ymd(dt); % year, month, day
            else
                if (isfield(these_data, 'yearcollected'))
                    year(d) = these_data.yearcollected;
                end
                month_direct = [];
                if (isfield(these_data, 'month'))
                    month_direct = str2double(these_data.month);
                end
                % eventDate comes in many different formats
                if (isfield(these_data, 'eventDate'))
                    [this_year, this_month, format_counter] = ...
                        parseDate(these_data.eventDate, format_counter);
                    if (isempty(this_month))
                        if (~isempty(month_direct))
                            month(d) = month_direct;
                        end
                    else                        
                        if (~isempty(month_direct) && ...
                                this_month ~= month_direct)
                            warning('mismatching months')
                            month_direct
                            this_month
                            month(d) = month_direct; % trust this one more
                        else
                            month(d) = this_month;
                        end
                    end
                    if (isfinite(this_year))
                        if (isfinite(year(d)) && this_year ~= year(d))
                            warning('different years: %s vs %d', ...
                                this_date, year(d))
                        end
                        year(d) = this_year; % always use this for consistency
                    elseif (length(these_data.eventDate) < 3 || ...
                            ~strncmp(these_data.eventDate,'---',3))
                        warning(['can''t parse eventDate: ', ...
                            these_data.eventDate])   
                        % keyboard
                    end
                end
            end
            for i = 1:length(anc_fields_string)
                if (isfield(these_data, anc_fields_string{i}))
                    eval(['tmp = these_data.', anc_fields_string{i}, ';'])
                    if (~isnumeric(tmp))
                        eval([anc_fields_string{i}, '(d) = nan;'])
                    else
                        eval([anc_fields_string{i}, ...
                            '(d) = str2num(these_data.', ...
                            anc_fields_string{i}, ');'])
                    end
                end
            end % i
            for i = 1:length(anc_fields_num)
                if (isfield(these_data, anc_fields_num{i}))
                    eval(['value = these_data.',anc_fields_num{i}, ';'])
                    if (~isempty(value))
                        eval([anc_fields_num{i}, '(d) = these_data.', ...
                            anc_fields_num{i}, ';']);
                    end
                end
            end
            for i = 1:length(anc_fields_char)
                if (isfield(these_data, anc_fields_char{i}))
                    eval(['value = these_data.',anc_fields_char{i}, ';'])
                    if (~isempty(value))
                        eval([anc_fields_char{i}, '{d} = these_data.', ...
                            anc_fields_char{i}, ';']);
                    end
                end
            end
            for i = 1:length(anc_string)
                eval(['no_val = isempty(', anc_string{i}, ');'])
                if (no_val && isfield(these_data, anc_string{i}))
                    eval([anc_string{i}, ' = these_data.', ...
                        anc_string{i}, ';'])
                end
            end % i
            for i = 1:length(anc_num)
                eval(['no_val = isnan(', anc_num{i}, ');'])
                if (no_val && isfield(these_data, anc_num{i}))
                    eval([anc_num{i}, ' = these_data.', ...
                        anc_num{i}, ';'])
                end
            end % i 
            if (isnan(minimumDepthInMeters(d)))
                depth(d) = maximumDepthInMeters(d);
            elseif (isnan(maximumDepthInMeters(d)))
                depth(d) = minimumDepthInMeters(d);
            else
                depth(d) = 0.5 * (minimumDepthInMeters(d) + ...
                    maximumDepthInMeters(d));
            end
        end % d
        % redefine longitude to match WOA convention (0..360)
        lon(lon<0) = lon(lon<0) + 360.;
        % this is a weird quirk in some data files: lon/lat/depth is all 0
        % mark them as bad data (real points are unlikely to have lon and
        % lat exactly 0)
        bad_idx = find(lon==0 & lat==0 & depth ==0);
        lon(bad_idx) = nan; 
        % fill in negative depth values for all missing depths, so that
        % depth=0 will be used for matching WOA values; but it's still
        % clear that they are actually missing
        no_depth_idx = find(isnan(depth));
        depth(isnan(depth)) = -10;
        % determine WOA temperature, salt, and oxygen for all values
        % that have all grid plus month values
        good_idx = find(~isnan(lon) & ~isnan(lat) & ~isnan(depth));
        if (~isempty(good_idx))
            longgrid = lon(good_idx);
            latgrid = lat(good_idx);
            depthgrid = depth(good_idx);
            yeargrid = year(good_idx);
            monthgrid = month(good_idx);
            mindepthgrid = minimumDepthInMeters(good_idx);
            maxdepthgrid = maximumDepthInMeters(good_idx);
            grid.xt2 = repmat(woa_grid.xt, [numel(longgrid) 1]);
            grid.yt2 = repmat(woa_grid.yt, [numel(longgrid) 1]);
            grid.zt2 = repmat(woa_grid.zt, [numel(longgrid) 1]);
            [~,minX] = min(abs(longgrid - grid.xt2),[],2);
            [~,minY] = min(abs(latgrid - grid.yt2),[],2);
            [~,minZ] = min(abs(depthgrid - grid.zt2),[],2);
            [mini,minZmin] = nanmin(abs(mindepthgrid - grid.zt2),[],2);
            minZmin(isnan(mini)) = nan;
            [mini,minZmax] = nanmin(abs(maxdepthgrid - grid.zt2),[],2);
            minZmax(isnan(mini)) = nan;
            temp = nan * ones(length(good_idx),12);
            salt = nan * ones(length(good_idx),12);
            o2 = nan * ones(length(good_idx),12);

            for i = 1:length(good_idx)
                % this is the index for the January value
                [idx, this_min_y, this_min_x] = find_idx_woa(WOA_mon.temp, ...
                    minZ(i), minY(i), minX(i), ...
                    lon(good_idx(i)), lat(good_idx(i)), woa_grid);
                minY(i) = this_min_y;
                minX(i) = this_min_x;
                temp(i,:) = WOA_mon.temp(idx:idx+11); % full year
                salt(i,:) = WOA_mon.salt(idx:idx+11);
                o2(i,:) = WOA_mon.o2(idx:idx+11);
            end % i
            
            for i = 1:length(anc_fields_multi)
                eval(['grid_',anc_fields_multi{i}, ' = ', ...
                    anc_fields_multi{i}, '(good_idx);'])
            end % i
            
            all_lon = [all_lon; longgrid];
            all_lat = [all_lat; latgrid];
            all_depth = [all_depth; depthgrid];
            all_year = [all_year; yeargrid];
            all_month = [all_month; monthgrid];
            all_temp = [all_temp; temp];
            all_salt = [all_salt; salt];
            all_o2 = [all_o2; o2];
            all_z = [all_z; minZ];
            all_y = [all_y; minY];
            all_x = [all_x; minX];
            all_zmin = [all_zmin; minZmin];
            all_zmax = [all_zmax; minZmax];
            for i = 1:length(anc_fields_string)
                eval(['all_',anc_fields_string{i}, ' = [ all_', ...
                    anc_fields_string{i}, '; grid_', ...
                    anc_fields_string{i}, '];']);
            end % i
            for i = 1:length(anc_fields_num)
                eval(['all_',anc_fields_num{i}, ' = [ all_', ...
                    anc_fields_num{i}, '; grid_', ...
                    anc_fields_num{i}, '];']);
            end % i
            for i = 1:length(anc_fields_char)
                eval(['empty_idx = find(cellfun(''isempty'', grid_', ...
                    anc_fields_char{i}, '));'])
                for e = 1:length(empty_idx)
                    eval(['grid_', anc_fields_char{i}, ...
                        '{empty_idx(e)} = ''UNKNOWN'';']);
                end % e
            end % i
            % for the first assignment, just copy from "grid_" to 
            % "all_
            if (is_first)
                for i = 1:length(anc_fields_char)
                    eval(['all_',anc_fields_char{i}, ' = grid_', ...
                        anc_fields_char{i}, ';']); % cell array
                end % i
                is_first = 0;
            else % subsequent files for same species
                for i = 1:length(anc_fields_char)
                    eval(['all_',anc_fields_char{i}, ' = cat(2, all_', ...
                        anc_fields_char{i}, ', grid_', ...
                        anc_fields_char{i}, '{:});']); % cell array
                end % i
            end
            clear temp salt o2;
        end    
        index = index + 1;
        fname = strrep(fname, num2str(index-1), num2str(index));        
    end
    
    if (~isempty(all_lon)) % don't create empty files
        % rename variables for output
        lon = all_lon; lat = all_lat; depth = all_depth;
        year = all_year; month = all_month;
        temp = all_temp; salt = all_salt; o2 = all_o2;
        z = all_z;
        y = all_y;
        x = all_x;
        zmin = all_zmin;
        zmax = all_zmax;
        % determine bathymetry from etopo
        bathymetry = nan * ones(size(all_x));
        for d = 1:length(all_x)
            [minx, idx] = min(abs(lon(d) - etopo.lon));
            [miny, idy] = min(abs(lat(d) - etopo.lat));
            if (minx > 0.01 || miny > 0.01)
                % some hard-coding here: grid resolution is 1 min, so nearest
                % point should never be more than 1/120 dg away, round up a bit
                warning('unexpected')
                keyboard
            end
            bathymetry(d) = etopo.topo(idy,idx);
        end % d        
        
        for i = 1:length(anc_fields_multi)
            eval([anc_fields_multi{i}, ' = all_', ...
                anc_fields_multi{i}, ';']);
        end % i
        %        
        fn_mat = [out_dir, strrep(files(f).name, '_1.json', '_mapped.mat')];
        save(fn_mat, 'lon', 'lat', 'depth', 'year', 'month', 'temp', ...
            'salt', 'o2', 'x', 'y', 'z', 'zmin', 'zmax', 'bathymetry');
        for i = 1:length(anc_fields_multi)
            save(fn_mat, anc_fields_multi{i}, '-append');
        end % i
        for i = 1:length(anc_fields_single)
            save(fn_mat, anc_fields_single{i}, '-append');
        end % i
    end
    clear temp salt o2 grid;
end % f
close(hw)

%% find nearest ocean cell if the original cell (input min_z/y/x) has
% nan values (=land) in a 1-cell circle around the original cell
% take the ocean cell that is nearest to the original cell if there are
%  multiple ocean cells within 1 cell from the original cell
function [idx, min_y, min_x] = find_idx_woa(woa_temp, min_z, min_y, min_x, ...
    this_lon, this_lat, woa_grid)
    idx = sub2ind(size(woa_temp), 1, min_z, ...
        min_y, min_x);
    temp = woa_temp(idx); % test with January, mask is same for all months
    if (isnan(temp))
        yrange = [max(min_y-1,1):min(min_y+1,180)];
        % FIXME not doing wraparounds for now
        xrange = [max(min_x-1,1):min(min_x+1,360)];
        temp_circle = squeeze(woa_temp(1,min_z,yrange,xrange));
        is_good_temp = isfinite(temp_circle);
        good_idx = find(is_good_temp);
        if (isempty(good_idx))
            %disp('no temp data in range')
            return
        end
        [x2d,y2d] = meshgrid(xrange,yrange);
        woa_lon = woa_grid.xt(xrange);
        woa_lat = woa_grid.yt(yrange);
        [lon2d,lat2d] = meshgrid(woa_lon,woa_lat);
        good_lon = lon2d(good_idx);
        good_lat = lat2d(good_idx);
        all_lon = this_lon * ones(length(good_idx)*2,1);
        all_lon(2:2:end) = good_lon;
        all_lat = this_lat * ones(length(good_idx)*2,1);
        all_lat(2:2:end) = good_lat;
        full_range = m_lldist(all_lon,all_lat); % includes duplicates
        all_range = full_range(1:2:end);
        [~,mini] = min(all_range);
        good_xx = x2d(good_idx);
        good_yy = y2d(good_idx);
        min_y = good_yy(mini);
        min_x = good_xx(mini);       
        idx = sub2ind(size(woa_temp), 1, min_z, ...
            min_y, min_x);
    end    
end
    
function [year, month, format_counter] = parseDate(date, format_counter)
% D/M/Y, M/D/Y, Y/M/D are the commonly used formats that have 3 numbers
% but there are many others (just year, year/month, year1-year2), extra
% strings ("circa", time stamps etc.)
year = [];
month = [];
% date string must contain numbers (otherwise, e.g.: "---")
if (~any(date >= '0' & date <= '9'))
    return;
end
months = {'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug'; ...
    'Sep';'Oct';'Nov';'Dec'};
% try to set up a standardized format
this_date = strrep(strrep(date, 'ca.', ''), 'circa', '');
this_date = strrep(this_date,'approx.','');
this_date = strrep(strrep(strrep(this_date, '/', '-'), ' ','-'),'\t','-');
this_date = strrep(strrep(this_date, '(', ''), ')', '');
this_date = strrep(strrep(this_date, '?', ''), '.', '-');
while (strcmp(this_date(1),'-'))
    this_date = this_date(2:end);
end
this_date = strtok(this_date,';');
this_date = strtok(this_date,'+');
if (isempty(findstr(this_date,'OCT')))
    this_date = strtok(this_date,'T');
end
while (strcmp(this_date(end),'-'))
    this_date = this_date(1:end-1);
end
%
[part1, remain1] = strtok(this_date, '-');
num1 = str2double(part1);
if (isnan(num1))
    % this happens if the first string doesn't represent a number, e.g.
    % if it is a string representation of a month ("Mar-91")
    if (all(isletter(part1)))
        for m = 1:12
            if (strncmpi(part1, months{m},3))
                month = m;
                [part2, ~] = strtok(remain1, '-');
                yr = str2double(part2);
                if (isfinite(yr))
                    year = yr;
                else
                    date
                    warning('CASE1 unexpected format')
                end
            end
        end % m
    else
        date
        warning('CASE2 unexpected format')
    end
elseif (isempty(remain1) || strncmp(remain1, '+00', 3) || ...
    (~any(remain1 >= '0' & remain1 <= '9') ) )
    % only year given
    year = num1;
else    
    [part2, remain2] = strtok(remain1, '-');
    % deal with something like this: '06T23:57:00+00:00', but 
    % don't modify a 3-letter month, e.g., '(18 AUG 1966)'
    if (length(part2) > 3 && isnan(str2double(part2(3))))
        part2 = part2(1:2);
    end
    num2 = str2double(part2);
    if (isnan(num2))
        % this could be day/month/year with month as string, e.g.
        % '(18 AUG 1966)'
        if (all(isletter(part2)))
            for m = 1:12
                if (strncmpi(part2, months{m},3))
                    month = m;
                    part3 = strtok(remain2, '-');
                    yr = str2double(part3);
                    if (isfinite(yr))
                        year = yr;
                    else
                        date
                        warning('CASE3 unexpected format')
                    end
                end
            end % m
        else
            date
            warning('CASE4 unexpected format')
        end
    elseif (isempty(remain2)) 
        % nothing left: only two numbers were given
        if (num2 >= 1 && num2 <= 12)
            % assume "year/month" format
            year = num1;
            month = num2;
        elseif (num2 >= num1)
            % assume "year1/year2" format, no month given
            date
            year = num1
            warning('two years assumed, using the first')
        elseif (num2 == 0)
            % bogus month, use only the year
            date
            year = num1
            warning('month 0 assumed')
        else
            date
            warning('CASE5 unexpected format')
        end
    else % date contains at least 2 numbers
        part3 = strtok(remain2, '-');
        num3 = str2double(part3);
        if (isnan(num3))
            % only two actual numbers, assume last part is garbage or the
            % day
            date
            warning('assuming year/month')
            year = num1;
            month = num2;
            %keyboard
        elseif (num1 < 1 || num2 < 1 || num3 < 1)
            if (num3 == 0) % this must be year 2000
                year = num3;
                month = num1;
                date
                warning('assuming m/d/y')
            elseif (num1 == 0) % this must be year 2000
                year = num1;
                month = num2;
                date
                warning('assuming y/m/d')
            else
                date
                warning('CASE6 unexpected format')
            end
        end
        % 3 actual numbers - use them only if it is quite clear which
        % format is used
        if (num1 > 31)
            year = num1;
            if (num3 > 12)
                % year/month/day
                month = num2;
            elseif (num2 > 12)
                % year/day/month
                month = num3;
            else % assume year/month/day (year/day/month would be
                % very unusual)
                month = num2;
            end
        else % num1 is in 1..31 range, could be year, month, or day
            if (num3 > 31)
                % that must be the year
                year = num3;
                if (num1 == num2)
                    month = num1; % same day and month
                elseif (num1 > 12)       
                    % day/month/year
                    month = num2;
                    format_counter(1) = format_counter(1) + 1;
                elseif (num2 > 12)
                    % month/day/year
                    month = num1;
                    format_counter(2) = format_counter(2) + 1;
                elseif (format_counter(1) > 10 && format_counter(2) == 0)
                    % assume day/month/year as it appears standard for
                    % this species
                    month = num2;
                elseif (format_counter(2) > 10 && format_counter(1) == 0)
                    % assume month/day/year as it appears standard for
                    % this species
                    month = num1;
                else
                    date
                    warning('CASE_A ambiguous format')
                end
            else % num1 and num3 are both in 1..31 range
                if (num2 > 31)
                    % year in the middle is unexpected
                    date
                    warning('CASE7 unexpected format')
                elseif (num2 > 12 && num3 > 12)
                    % month must be first number, so year the last
                    month = num1;
                    year = num3;
                else % e.g., 10/4/05
                     date
                     warning('CASE_B ambiguous format')
                end
            end
        end
    end
end

% deal with non-Y2K compliant data
if (year < 1000)
    % assume that small numbers represent the 21st century, others the 20th
    if (year < 22) % need to update this every year (currently set for 2021)
        year = year + 2000;
    else
        year = year + 1900;
    end
end

if (~isempty(month) && (month < 1 || month > 12))
    % non-parsable and bad data
    fprintf('Bad date (%s), month: %d\n', date, month)
    month = []
end


end
