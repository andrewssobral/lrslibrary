function GCO_UnitTest
% GCO_UnitTest   Compile, load, and test the GCO_MATLAB library.
%    GCO_UnitTest will make sure the wrapper compiles on the target 
%    platform and then exercises the library to look for silly bugs.

    function Assert(cond,msg)   % for older MATLAB without assert()
        if (exist('assert') == 5)
            if (nargin < 2)
                assert(cond);
            else
                assert(cond,msg);
            end
        elseif (~cond)
            if (nargin < 2)
                msg = 'Assertion failed';
            end
            error(msg);
        end
    end

    function DoSetDataCost(h,dc,iter)
        if (iter == 1)
            % Set data costs as dense matrix
            GCO_SetDataCost(h,dc);
        else
            % Use the sparse mechanism to set dense data costs,
            % to verify that tests all give the same result.
            % Note that this is not a good test when the number
            % of sites is small, since only some of the sparse datacost 
            % code path will get exercised.
            for label=1:size(dc,1)
                ids = find(dc(label,:) < 100000);
                GCO_SetDataCost(h,[ids; dc(label,ids)],label);
            end
        end
    end
    
sc = [0 1 2 3 3 3 3 3 3;
      1 0 1 2 3 3 3 3 3;
      2 1 0 1 2 3 3 3 3;
      3 2 1 0 1 2 3 3 3;
      3 3 2 1 0 1 2 3 3;
      3 3 3 2 1 0 1 2 3;
      3 3 3 3 2 1 0 1 2;
      3 3 3 3 3 2 1 0 1;
      3 3 3 3 3 3 2 1 0;];  % truncated linear

GCO_BuildLib; disp('BuildLib PASSED');
GCO_LoadLib;  disp('LoadLib PASSED');

% Basic tests with no Create/Delete
caught=false; try GCO_Delete(10);          catch, caught=true; end, Assert(caught,'Expected an exception');
caught=false; try h = GCO_Create(1,1);     catch, caught=true; end, Assert(caught,'Expected an exception');
caught=false; try h = GCO_Create(0,2);     catch, caught=true; end, Assert(caught,'Expected an exception');
h1 = GCO_Create(4,8);
h2 = GCO_Create(10,5);
Assert(all(GCO_ListHandles == [h1; h2]));
GCO_Delete(h1);
caught=false; try GCO_ComputeEnergy(h1);   catch, caught=true; end, Assert(caught,'Expected an exception');
caught=false; try GCO_Delete(h1);          catch, caught=true; end, Assert(caught,'Expected an exception');
Assert(all(GCO_ListHandles == [h2]));
Assert(length(GCO_GetLabeling(h2)) == 10);
GCO_Delete(h2);
Assert(isempty(GCO_ListHandles));
caught=false; try GCO_Delete(h2);          catch, caught=true; end, Assert(caught,'Expected an exception');
disp('Create/Delete PASSED');
disp('ListHandles PASSED');

% Test Get/SetLabel when no optimization is used
h = GCO_Create(4,3);
Assert(GCO_ComputeEnergy(h) == 0);
l = (mod((1:4),3)+1)';
caught=false; try GCO_SetLabeling(h,l(1:end-1)); catch, caught=true; end, Assert(caught,'Expected an exception');
caught=false; try GCO_SetLabeling(h,[l 1]);      catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_SetLabeling(h,l);
Assert(all(GCO_GetLabeling(h) == l));
Assert(all(GCO_GetLabeling(h,3) == l(3)));
Assert(all(GCO_GetLabeling(h,2,2) == l(2:3)));
caught=false; try GCO_GetLabeling(h,0,1);       catch, caught=true; end, Assert(caught,'Expected an exception');
caught=false; try GCO_GetLabeling(h,1,5);       catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_SetLabeling(h,4-l);
Assert(all(GCO_GetLabeling(h) == 4-l));
GCO_Delete(h);
disp('Get/SetLabeling PASSED');

% Test with NO costs 
h = GCO_Create(4,3);
Assert(GCO_ComputeEnergy(h) == 0);
GCO_Expansion(h);
Assert(GCO_ComputeEnergy(h) == 0);
GCO_SetLabeling(h,4-GCO_GetLabeling(h));
Assert(GCO_ComputeEnergy(h) == 0);
GCO_Expansion(h);
GCO_Delete(h);
disp('Expansion-000 PASSED');

for iter=1:2
if (iter == 2)
    fprintf('SetDataCost-Sparse...'); 
    h = GCO_Create(3,2);
    GCO_SetDataCost(h,[1 1 1; 2 2 2]); % set dense, then make sure sparse is not allowed afterwards
    caught=false; try GCO_SetDataCost(h,[1 2 3; 5 5 5],1);  catch, caught=true; end, Assert(caught,'Expected an exception');
    GCO_Delete(h);
%     h = GCO_Create(3,2);
%     GCO_SetDataCost(h,[1 2 3; 5 6 7],1); % set sparse, then make sure dense is not allowed afterwards
%     GCO_SetDataCost(h,[1; 100 ],2); % make sure label cost can be replaced
%     Assert(GCO_ComputeEnergy(h) == 18);
%     GCO_SetDataCost(h,[1 2 3; 10 11 12],1); % make sure label cost can be replaced
%     Assert(GCO_ComputeEnergy(h) == 33);
%     % make sure unsorted ids will fail order is checked
%     caught=false; try GCO_SetDataCost(h,[1 1 2; 9 9 9],1);  catch, caught=true; end, Assert(caught,'Expected an exception');
%     % make sure dense can be used after sparse (though not vice versa)
%     GCO_SetDataCost(h,[1 1 1; 2 2 2]);
%     GCO_Delete(h);
end

% Test Expansion with DATA cost only
h = GCO_Create(4,9);
dc = [1 2 5 8 4 2 3 7 9; 
      3 1 2 5 4 5 5 5 5;
      5 5 5 5 4 5 2 1 3;
      9 7 3 2 4 8 5 2 1;]';
if (iter == 1)
    caught=false; try GCO_SetDataCost(h,[dc [0 0 0 0]']);    catch, caught=true; end, Assert(caught,'Expected an exception');
    caught=false; try GCO_SetDataCost(h,dc(:,1:end-1));      catch, caught=true; end, Assert(caught,'Expected an exception');
    caught=false; try GCO_SetDataCost(h,dc(1:end-1,:));      catch, caught=true; end, Assert(caught,'Expected an exception');
else
    caught=false; try GCO_SetDataCost(h,[1:4; dc(1,:)],0);   catch, caught=true; end, Assert(caught,'Expected an exception');
    caught=false; try GCO_SetDataCost(h,[1:4; dc(1,:)],10);  catch, caught=true; end, Assert(caught,'Expected an exception');
    caught=false; try GCO_SetDataCost(h,[0:3; dc(1,:)],1);   catch, caught=true; end, Assert(caught,'Expected an exception');
    caught=false; try GCO_SetDataCost(h,[2:5; dc(1,:)],1);   catch, caught=true; end, Assert(caught,'Expected an exception');
    caught=false; try GCO_SetDataCost(h,[1 2 2;dc(1,1:3)],1);catch, caught=true; end, Assert(caught,'Expected an exception');
end
DoSetDataCost(h,dc,iter);
GCO_SetLabeling(h,[3 3 3 3]);
Assert(GCO_ComputeEnergy(h) == sum(dc(3,:)));
GCO_SetLabeling(h,[1 2 3 4]);
Assert(GCO_ComputeEnergy(h) == dc(1,1)+dc(2,2)+dc(3,3)+dc(4,4));
GCO_SetLabeling(h,[5 5 5 5]);
GCO_ExpandOnAlpha(h,8);
Assert(all(GCO_GetLabeling(h) == [5 5 8 8]'));
GCO_ExpandOnAlpha(h,3);
Assert(all(GCO_GetLabeling(h) == [5 3 8 8]'));
GCO_Expansion(h);
Assert(all(GCO_GetLabeling(h) == [1 2 8 9]'));
Assert(GCO_ComputeEnergy(h) == dc(1,1)+dc(2,2)+dc(8,3)+dc(9,4));
GCO_Delete(h);
if (iter==1), disp('Expansion-D00 PASSED'); end

% Test DATA+SMOOTH cost 
h = GCO_Create(4,9);
dc = [1 2 5 8 4 2 3 7 9; 
      3 1 1 5 4 5 5 5 5;
      5 5 5 5 4 5 1 1 3;
      9 7 3 2 4 8 5 2 1;]';
DoSetDataCost(h,dc,iter);
caught=false; try GCO_SetSmoothCost(h,sc(:,1:end-1));        catch, caught=true; end, Assert(caught,'Expected an exception');
caught=false; try GCO_SetSmoothCost(h,sc(1:end-1,:));        catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_SetSmoothCost(h,sc);
caught=false; try GCO_SetNeighbors(h,eye(4));              catch, caught=true; end, Assert(caught,'Expected an exception');
caught=false; try GCO_SetNeighbors(h,zeros(5));            catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_SetNeighbors(h,[0 2 0 0;
                    0 0 1 0;
                    0 0 0 2;
                    0 0 0 0]);
GCO_SetLabeling(h,[3 3 3 3]);
Assert(GCO_ComputeEnergy(h) == sum(dc(3,:)));
GCO_SetLabeling(h,[1 2 4 5]);
Assert(GCO_ComputeEnergy(h) == dc(1,1)+dc(2,2)+dc(4,3)+dc(5,4) + 6);
GCO_SetLabeling(h,[5 5 5 5]);
GCO_Expansion(h);
Assert(all(GCO_GetLabeling(h) == [2 2 8 8]'));
GCO_Delete(h);
if (iter==1) disp('Expansion-DS0 PASSED'); end

% Test DATA+LABEL cost 
h = GCO_Create(4,9);
dc = [1 2 5 8 4 2 3 7 9; 
      3 1 3 5 4 3 5 5 5;
      5 5 5 5 4 5 1 1 3;
      9 7 3 2 4 8 5 2 1;]';
lc = [9 9 1 1 1 2 1 9 9];
DoSetDataCost(h,dc,iter);
caught=false; try GCO_SetLabelCost(h,[lc 10]);            catch, caught=true; end, Assert(caught,'Expected an exception');
caught=false; try GCO_SetLabelCost(h,lc(1:end-1));        catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_SetLabelCost(h,lc);
GCO_SetLabeling(h,[3 3 3 3]);
Assert(GCO_ComputeEnergy(h) == sum(dc(3,:))+lc(3));
GCO_SetLabeling(h,[1 2 4 5]);
Assert(GCO_ComputeEnergy(h) == dc(1,1)+dc(2,2)+dc(4,3)+dc(5,4) + sum(lc([1 2 4 5])));
GCO_SetLabeling(h,[5 5 5 5]);
GCO_Expansion(h);
Assert(all(GCO_GetLabeling(h) == [7 3 7 3]'));
GCO_Delete(h);

% Test when NumLabels < NumSites and make sure greedy doesn't crash when all labels 
% are added. This test thanks to Yangyan Li.
h = GCO_Create(4,3);
DoSetDataCost(h,[1 4 9; 5 2 5; 6 1 3; 5 7 1;]',iter);
GCO_SetLabelCost(h,[1 1 1]);
GCO_Expansion(h);
GCO_Delete(h);


% Now do the same test, except add label costs to subsets of labels, not
% just individual labels
h = GCO_Create(4,9);
DoSetDataCost(h,dc,iter);
GCO_SetLabelCost(h,lc);
GCO_SetLabelCost(h,3,[1 3 4 5 6]);
GCO_SetLabelCost(h,4,[4 8 9]);
GCO_SetLabeling(h,[3 3 3 3]);
Assert(GCO_ComputeEnergy(h) == sum(dc(3,:))+lc(3)+3);
GCO_SetLabeling(h,[5 5 8 9]);
Assert(GCO_ComputeEnergy(h) == dc(5,1)+dc(5,2)+dc(8,3)+dc(9,4) + sum(lc([5 8 9]))+3+4);
GCO_SetLabeling(h,[1 2 4 5]);
Assert(GCO_ComputeEnergy(h) == dc(1,1)+dc(2,2)+dc(4,3)+dc(5,4) + sum(lc([1 2 4 5]))+3+4);
GCO_SetLabeling(h,[5 5 5 5]);
GCO_Expansion(h);
Assert(all(GCO_GetLabeling(h) == [7 7 7 7]'));
GCO_SetLabelCost(h,3,[2]);
GCO_Expansion(h);
Assert(all(GCO_GetLabeling(h) == [2 2 7 7]'));
GCO_Delete(h);
if (iter==1), disp('Expansion-D0L PASSED'); end

% Test DATA+SMOOTH+LABEL cost 
h = GCO_Create(4,9);
dc = [1 2 5 8 4 2 3 7 9; 
      3 1 3 5 4 2 5 5 5;
      5 5 5 5 5 5 1 1 3;
      9 7 3 2 4 8 5 2 1;]';
lc = [1 9 1 1 1 1 1 9 9];
GCO_SetSmoothCost(h,sc);
GCO_SetNeighbors(h,[0 2 0 0;
                    0 0 1 0;
                    0 0 0 2;
                    0 0 0 0]);
DoSetDataCost(h,dc,iter);
GCO_SetLabelCost(h,lc);
GCO_SetLabeling(h,[3 3 3 3]);
Assert(GCO_ComputeEnergy(h) == sum(dc(3,:))+lc(3));
GCO_SetLabeling(h,[1 2 4 5]);
Assert(GCO_ComputeEnergy(h) == dc(1,1)+dc(2,2)+dc(4,3)+dc(5,4) + sum(lc([1 2 4 5])) + 6);
GCO_SetLabeling(h,[5 5 5 5]);
GCO_Expansion(h);
Assert(all(GCO_GetLabeling(h) == [6 6 7 7]'));
GCO_Delete(h);
% Now do the same test, except add label costs to subsets of labels, not
% just individual labels
h = GCO_Create(4,9);
GCO_SetSmoothCost(h,sc);
GCO_SetNeighbors(h,[0 2 0 0;
                    0 0 1 0;
                    0 0 0 2;
                    0 0 0 0]);
DoSetDataCost(h,dc,iter);
GCO_SetLabelCost(h,lc);
GCO_SetLabelCost(h,3,[1 3 4 5 6]);
GCO_SetLabelCost(h,4,[5 8 9]);
GCO_SetLabeling(h,[3 3 3 3]);
Assert(GCO_ComputeEnergy(h) == sum(dc(3,:))+lc(3)+3);
GCO_SetLabeling(h,[5 5 8 9]);
Assert(GCO_ComputeEnergy(h) == dc(5,1)+dc(5,2)+dc(8,3)+dc(9,4) + sum(lc([5 8 9]))+3+4 + 5);
GCO_SetLabeling(h,[1 2 4 5]);
Assert(GCO_ComputeEnergy(h) == dc(1,1)+dc(2,2)+dc(4,3)+dc(5,4) + sum(lc([1 2 4 5]))+3+4 + 6);
GCO_SetLabeling(h,[5 5 5 5]);
GCO_Expansion(h);
Assert(all(GCO_GetLabeling(h) == [7 7 7 7]'));
GCO_SetLabelCost(h,0,[7]);   % Try replacing an existing labelcost
GCO_SetLabeling(h,[3 3 3 3]);
GCO_Expansion(h);
Assert(all(GCO_GetLabeling(h) == [7 7 7 7]'));
[E D S L] = GCO_ComputeEnergy(h);
Assert(L == 0);
GCO_SetLabelCost(h,4,[7]);
GCO_SetLabelCost(h,1,[1 3 4 5 6]);  % Try replacing a labelcost subset
GCO_Expansion(h);
Assert(all(GCO_GetLabeling(h) == [6 6 4 4]'));
[E D S L] = GCO_ComputeEnergy(h);
Assert(L == 3);
GCO_SetSmoothCost(h,sc*3);   % Try replacing smoothcost
GCO_SetLabeling(h,[5 5 7 9]);
Assert(GCO_ComputeEnergy(h) == dc(5,1)+dc(5,2)+dc(7,3)+dc(9,4) + sum(lc([5 9]))+1+4+4 + 6*3);
GCO_Delete(h);
if (iter==1), disp('Expansion-DSL PASSED'); end

% Test NON-METRIC SMOOTH cost, and make sure Expansion raises an exception
% so that users do not accidentally get meaningless results.
h = GCO_Create(4,9);
dc = [1 1 1 1 1 1 1 1 1; 
      1 1 1 1 1 1 1 1 1;
      1 1 1 1 1 1 1 1 1;
      1 1 1 1 1 1 1 1 1;]';
DoSetDataCost(h,dc,iter);
sc_nonmetric = sc.*sc;  % truncated quadratic
GCO_SetSmoothCost(h,sc_nonmetric);
GCO_SetNeighbors(h,[0 2 0 0;
                    0 0 1 0;
                    0 0 0 2;
                    0 0 0 0]);
GCO_SetLabeling(h,[5 5 7 4]);
caught=false; try GCO_ExpandOnAlpha(h,5);           catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_Delete(h);
if (iter==1), disp('Expansion-NonMetricWarning PASSED'); end


if iter==1
% Test Swap 
h = GCO_Create(4,9);
dc = [1 2 5 8 4 2 3 7 9; 
      3 1 1 5 4 5 5 5 5;
      5 5 5 5 4 5 1 1 3;
      9 7 3 2 4 8 5 2 1;]';
GCO_SetDataCost(h,dc);
GCO_SetLabeling(h,[5 5 5 5]);
GCO_Swap(h);
Assert(all(GCO_GetLabeling(h) == [1 2 7 9]'));
GCO_Delete(h);
disp('Swap-D0 PASSED');
h = GCO_Create(4,9);
GCO_SetDataCost(h,dc);
GCO_SetSmoothCost(h,sc);
GCO_SetNeighbors(h,[0 2 0 0;
                    0 0 1 0;
                    0 0 0 2;
                    0 0 0 0]);
GCO_SetLabeling(h,[5 5 5 5]);
GCO_Swap(h);
Assert(all(GCO_GetLabeling(h) == [2 2 8 8]'));
GCO_Delete(h);
disp('Swap-DS PASSED');
end

% Make sure to generate errors if dangerously large coefficients are used
dc = [1 1 1 1 1 1 1 1 1; 
      1 1 1 1 1 1 1 1 1;
      1 1 1 1 1 1 1 1 10000001; % huge capacity, relative to int32
      1 1 1 1 1 1 1 1 1;]';
lc = [1 1 0 0 0 0 0 0 10000001];
h = GCO_Create(4,9);
caught=false; try GCO_SetLabelCost(h,lc);           catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_SetDataCost(h,dc);   % test BadLabelCost, and then BadDataCost code path
caught=false; try GCO_ExpandOnAlpha(h,9);           catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_Delete(h);

h = GCO_Create(4,9);
GCO_SetDataCost(h,dc);
GCO_SetSmoothCost(h,sc); % test BadDataCost+GoodSmoothCost code path
caught=false; try GCO_ExpandOnAlpha(h,9);           catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_Delete(h);

h = GCO_Create(4,9);
DoSetDataCost(h,min(dc,1),iter);
GCO_SetSmoothCost(h,sc.*10); % test GoodDataCost+BadSmoothCost code path
GCO_SetNeighbors(h,[0 2 0 0;
                    0 0 1 0;
                    0 0 0 10000001;
                    0 0 0 0]);
GCO_SetLabeling(h,[1 1 1 1]);
caught=false; try GCO_ExpandOnAlpha(h,5);           catch, caught=true; end, Assert(caught,'Expected an exception');
GCO_Delete(h);
if (iter==1), disp('IntegerOverflowWarnings PASSED'); end
if (iter==2), fprintf('PASSED\n'); end


% Test MEDIUM SCALE problems and make sure dense/sparse data costs
% result in the exact same solution
rand('twister', 987); % get the same random stream each time
wd = 64; ht = 48;
dc = int32(rand([500,wd*ht])*1000);
dc(dc > 100) = 100000; % prune out about 10% of possible labels
dc(1,:) = 1000; % but be sure to allow at least one label per variable
nb = sparse(wd*ht,wd*ht);
for y=1:ht % set up a grid-like neighbourhood, arbitrarily
    for x=1:wd
        if (x < wd), nb((y-1)*wd+x,(y-1)*wd+x+1) = 30; end
        if (y < ht), nb((y-1)*wd+x, y   *wd+x  ) = 30; end
    end
end

% Test greedy
if (iter==1), fprintf('MediumScale-D0L...'); end
if (iter==2), fprintf('MediumScale-D0L-Sparse...'); end
h = GCO_Create(wd*ht,size(dc,1));
DoSetDataCost(h,dc,iter);
GCO_SetLabelCost(h,3000);
tic; GCO_Expansion(h); greedytime = toc;
[E D S L] = GCO_ComputeEnergy(h);
if (exist('EDLdense'))
    Assert(GCO_ComputeEnergy(h) == EDLdense);  % Test D+L costs
else
    EDLdense = GCO_ComputeEnergy(h); % remember for next time, to compare with "sparse data cost" solution
end
GCO_Delete(h);
if (iter==1), fprintf('PASSED        (%.3fsec greedy)\n',greedytime); end
if (iter==2), fprintf('PASSED (%.3fsec greedy    w/ 10%% of sites feasible)\n',greedytime); end


% Test expansion
if (iter==1), fprintf('MediumScale-DSL...'); end
if (iter==2), fprintf('MediumScale-DSL-Sparse...'); end
h = GCO_Create(wd*ht,size(dc,1));
DoSetDataCost(h,dc,iter);
GCO_SetLabelCost(h,1000);
GCO_SetLabelCost(h,5000,2:50);
GCO_SetLabelCost(h,10000,100:200);
GCO_SetNeighbors(h,nb);
tic; GCO_Expansion(h); exptime = toc;

[E D S L] = GCO_ComputeEnergy(h);
if (exist('EDSLdense'))
    Assert(GCO_ComputeEnergy(h) == EDSLdense);  % Test D+S+L costs
else
    EDSLdense = GCO_ComputeEnergy(h); % remember for next time, to compare with "sparse data cost" solution
end
GCO_Delete(h);
if (iter==1), fprintf('PASSED        (%.3fsec expansion)\n',exptime); end
if (iter==2), fprintf('PASSED (%.3fsec expansion w/ 10%% of sites feasible)\n',exptime); end
end

Assert(isempty(GCO_ListHandles)); % expect there to be no gc handles remaining

end
