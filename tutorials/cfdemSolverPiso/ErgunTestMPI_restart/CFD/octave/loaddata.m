function [data,colname]=loaddata(filename,columns,headerlines)

% data = loaddata(filename,columns,headerlines)
% data = loaddata(filename,0,headerlines): automatc detection of 
%           columns by word count in first headerline or first data line
%           (if headerline==0), separated by blanks or tabs
%           columns=-1: write colum assignment code lines
%              works only if headerline>0 and colnames separated by 
%              not more than 1 blank!
%           columns=-2: columns separated by 2 blanks (fluent/scheme auswertung-instat.scm)
% [data, colname] = loaddata(filename,0,headerlines) 
%           colname = cell string array of column names from headerline

fprintf(1,'loading %s ...  ', filename);
f=fopen(filename,'r');
if f==-1
  fprintf(1,'\n*** error: could not open "%s" ...\n', filename);
  data=[];
else
  for i=1:headerlines
    if i==1 s=fgets(f); else fgets(f); end
  end
  if headerlines==0
    s=fgets(f); frewind(f);
  end
  if columns<=0 % & headerlines>=1
    fprintf(1,'\n');
    fprintf(1,'  %s',s);
    pos = findstr(sprintf('\t'),s); % trennzeichenpositionen in zeile suchen
    if length(pos)>0
      endpos=length(s);
      %if isspace(s(endpos-1)) endpos=endpos-1; end
      pos = [0 pos(1,:) endpos]; % anfangs- und endposition hinzufügen
    elseif columns==-2 % spalteneinträge durch 2 blanks getrennt (fluent/scheme)
      pos=[0];
      for i=1:length(s)-1
        if isspace(s(i)) & isspace(s(i+1))
          pos(end+1)=i+1;
        end
      end
      pos(end+1)=length(s);      
    else % spalteneinträge durch blanks getrennt     
      %pos = findstr(' ',s); % trennzeichenpositionen in zeile suchen
      pos=[];
      word=0;
      for i=1:length(s)
        if isspace(s(i))
          word=0;
        elseif ~word
          word=1;
          pos(end+1)=i-1;
        end
      end
      pos(end+1)=length(s);
    end
    if headerlines>=1
      for i=1:length(pos)-1 % alle spalten
        colname{i}=s(pos(i)+1:pos(i+1)-1);
        if columns==-1 fprintf(' = data(strmatch(''%s'',colname),:); %% column %d\n', colname{i},i); end
        if columns==-2 fprintf('%% column %d: %s\n', i, colname{i}); end
      end   
    end
    columns = length(pos)-1;
    fprintf(1,'  total: %d columns ',columns);
  end
  data=fscanf(f,'%f',[columns,inf]);
  fclose(f);
  fprintf(1,'done.\n');
end
