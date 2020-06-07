
% %% Task1 %%
% %sudoku_propagate(+,-)
% sudoku_propagate(Instance, List):-
%     Instance = sudoku(SqrtN, Hints),
%     N is SqrtN**2,
%     length(Board,N),
%     initiate_board(Board,N),
%     put_hints(Hints,Board),
%     % List is Board.
%     check_rows_constraints(Board,1,List_rows,N),

%     check_cols_constraints(Board,List_cols,N),
%     check_boxes_constraints(Board,N,SqrtN,Boxes_list),

%     append(List_rows,List_cols,T_list),
%     append(T_list,Boxes_list,List).

% initiate_board([],_N).
% initiate_board([Row|Rows],N):-
%     length(Row,N),
%     initiate_board(Rows,N).

% % check_cols(N,_Board,N,_List).
% % check_cols(Row_count,Board,N,List),
% %     nth1(Row_count,Board,Row),

% check_rows_constraints([],_Row_index,[],_N).
% check_rows_constraints([Row|Rest],Row_index,Rows_list,N):-
%     Next_row is Row_index+1,
%     (check_row(Row,Row_index,N,Constraint)->
%     check_rows_constraints(Rest,Next_row,Rest_list,N),
%     append([Constraint],Rest_list,Rows_list);
%     check_rows_constraints(Rest,Next_row,Rows_list,N)).


% % check_row([],_Row_list,[],Col,N).
% % check_row([A|As],Row_list,Nums_in_row,Col,Vars,N):-
% %     % Ncol is Col+1,
% %     \+Vars>1,
% %     nonvar(A)-> Num=[A], N_vars is Vars; (Num=[], N_vars is Vars+1),
% %     check_row(As,)
% check_row(Row,Row_index,N,Constraint):-
%     % between(1,N,Val),
%     once(findnsols(2,Val,try_row(Row,Val,N),Sols)),
%     length(Sols,1),
%     Sols=[Sol],
%     put_in_row(Row,Row_index,Sol,1,Constraint).
%     % Out = Row.

% try_row(Row,Val,N):-
%     between(1,N,Val),
%     once(not_in_row(Row,Val)).
% not_in_row([],_Val).
% not_in_row([A|As],Val):-
%     (nonvar(A)-> A\=Val; A is Val),
%     not_in_row(As,Val).    
% put_in_row([A|As],Row_index,Val,Col_index,Constraint):-
%     (var(A)-> 
%     A is Val,
%     Constraint = (cell(Row_index,Col_index)=Val)
%     ;
%     Next_Col is Col_index+1,
%     put_in_row(As,Row_index,Val,Next_Col,Constraint)).

% % %%%
% % check_cols_constraints([],_Row_index,[],_N).
% check_cols_constraints(Rows,Cols_list,N):-
%     transpose(Rows,Cols),
%     check_rows_constraints(Cols,1,Cols_list,N),
%     transpose(Cols,Rows).

% check_boxes_constraints(Board,N,SqrtN,Boxes_list):-
%     blocks_to_rows(Board,N,SqrtN,Boxes),
%     check_rows_constraints(Boxes,1,Boxes_list,N).
%     % blocks_to_rows(Boxes,N,SqrtN,Board).

% user:file_search_path(sat, '/home/hadar/Desktop/plsatsolver_src/satsolver').
% :- use_module(sat(satsolver)).
%*************sudoku_propagate**************
% remove duplicates taken from https://stackoverflow.com/questions/39435709/how-to-remove-duplicates-from-a-list-in-swi-prolog
delete3(_,[],[]).
delete3(X,[X|T],R):- delete3(X,T,R).
delete3(X,[H|T],[H|R]) :- delete3(X,T,R).

remove_duplicates([],[]).
remove_duplicates([H|T], [H|R]) :- 
    member(H,T),!,
    delete3(H,T,R1),
    remove_duplicates(R1,R).

remove_duplicates([H|T],[H|R]):-
    remove_duplicates(T,R).

sudoku_propagate(Instance, List):-
    sudoku_propagate(Instance,[],List1,[cell(1,1)=_],0,_ExpC),
    append(List1,List).

% i check by the flag if in the previous rownd some caonstraints has added
sudoku_propagate(Instance,PrevL,List,[cell(0,0)=A],Flag,ExpC):-
Flag>0->sudoku_propagate(Instance,PrevL,List,[cell(1,1)=A],0,ExpC);
Flag=:=0->List=[],ExpC=[].
 
sudoku_propagate(sudoku(Sqr,Instance),PrevL,List,Pos,Flag,ExpC):-
    (Pos\= [cell(0,0)=_],
    append(Instance,PrevL,AllList),
    \+allready_checked(AllList,Pos),
    constraints(AllList,Sqr,Pos,Constraints_list,ExplainCons),!,
    remove_duplicates(Constraints_list,Cons_nodup),!,
    length(Cons_nodup,N_curr),
    pow(Sqr,2,N),
    N_1 is N-1,
    find_next_pos(AllList,Pos,NextPos,N),
    N_1==N_curr)->                                              %if ther is n-1 constraints
    (find_assignment(Cons_nodup,N,Pos),
    ExpC=[ExplainCons->Pos|ResExp],
    Nflag is Flag+1,
    append(PrevL,Pos,NewPrev),
    List=[Pos|ResList],
    sudoku_propagate(sudoku(Sqr,Instance),NewPrev,ResList,NextPos,Nflag,ResExp));
    (Pos\= [cell(0,0)=_],
    append(Instance,PrevL,AllList),
    pow(Sqr,2,N),
    find_next_pos(AllList,Pos,NextPos,N),
    sudoku_propagate(sudoku(Sqr,Instance),PrevL,List,NextPos,Flag,ExpC)).

%find next position to continue
find_next_pos(AllList,[cell(I1,J1)=A1],Pos,N):-
    J1=:=N,I1=:=N->Pos=[cell(0,0)=_B];          %if we are in the last cell, flag for finish (0,0)
    J1=:=N-> NextJ is 1 ,NextI is I1+1,(\+allready_checked(AllList,cell(NextI,NextJ)=_)-> Pos=[cell(NextI,NextJ)=_];
    (find_next_pos(AllList,[cell(NextI,NextJ)=A1],Pos,N)));   %if we in the last col we move to next row
    J1<N-> Next is J1+1 ,(\+allready_checked(AllList,[cell(I1,Next)=A1])-> Pos=[cell(I1,Next)=_C];
    (find_next_pos(AllList,[cell(I1,Next)=_],Pos,N))).
    

%find the assighnment
find_assignment(Cons_nodup,N,[cell(_I2,_J2)=A2]):-
    numlist(1,N,L),
    delete_elements(Cons_nodup,L,A2).

delete_elements([X|Cons_nodup],L,A2):-
    delete(L,X,Res),
    delete_elements(Cons_nodup,Res,A2).
delete_elements([],[A2],A2).    

%find the constraints for the current cell   
constraints([cell(I1,J1)=A1|Res1],Sqr,[cell(I2,J2)=A2],[A1|Constraints_list],[cell(I1,J1)=A1|ExplainCons]):-
    (I1==I2;J1==J2;find_box(I1,J1,I2,J2,Sqr)),             %check if it is in the row or col or box
    constraints(Res1,Sqr,[cell(I2,J2)=A2],Constraints_list,ExplainCons).

constraints([],_Sqr,_Pos,[],[]).
constraints([cell(I1,J1)=_A1|Res1],Sqr,[cell(I2,J2)=A2],Constraints_list,ExplainCons):-
    (I1\=I2,J1\=J2,\+find_box(I1,J1,I2,J2,Sqr)),             %check if it is in the row or col or box
    constraints(Res1,Sqr,[cell(I2,J2)=A2],Constraints_list,ExplainCons).

find_box(I1,J1,RowP,ColP,Sqr):-
    find_box_helper(1,Sqr,RowP,Sqr,Row_box),
    find_box_helper(1,Sqr,ColP,Sqr,Col_box),
    member(I1,Row_box),member(J1,Col_box).

find_box_helper(Start,End,I1,Sqr,Row_box):-
    (numlist(Start,End,L1),member(I1,L1))->Row_box=L1;
    NewStart is Start+Sqr,NewEnd is End+Sqr ,find_box_helper(NewStart,NewEnd,I1,Sqr,Row_box). 

%finde if we allready have assigment for the current cell
allready_checked([],[_Pos]):-false.
allready_checked([cell(I1,J1)=_A1|_Res1],[cell(I2,J2)=_A2]):-
    I1==I2,J1==J2.

allready_checked([cell(I1,J1)=_A1|Res1],[cell(I2,J2)=A2]):-
        (I1\=I2;J1\=J2) -> allready_checked(Res1,[cell(I2,J2)=A2]).
    
%*****************sudoku_propagate_explain(Instance, Explain)***********************
sudoku_propagate_explain(sudoku(Sqr, Hints), Explain):-
    sudoku_propagate(sudoku(Sqr, Hints),[],_List1,[cell(1,1)=_A],0,Exp),
    pow(Sqr,2,N),
    no_dup_explain(Exp,Explain,N).

%remove duplcates from the explain
no_dup_explain([],[],_N).
no_dup_explain([Cons_list->Pos|Res1],No_dup_exp,N):-
    length(Cons_list,Len),
    Len\=N-1,numlist(1,N,NumList),no_dup(Cons_list,New_ConsList,NumList),
    No_dup_exp=[New_ConsList->Pos|Res2],no_dup_explain(Res1,Res2,N). 

no_dup_explain([Cons_list->Pos|Res1],No_dup_exp,N):-
    length(Cons_list,Len),
    Len==N-1,No_dup_exp=[Cons_list->Pos|Res2],no_dup_explain(Res1,Res2,N).


no_dup([],[],_).
no_dup([cell(P1,P2)=A|Cons_list],New_ConsList,NumList):-
    (member(A,NumList)->
    delete(NumList,A,NewNumlist),
    New_ConsList =[cell(P1,P2)=A|Res],
    no_dup(Cons_list,Res,NewNumlist));
    \+member(A,NumList)->no_dup(Cons_list,New_ConsList,NumList).    

% the code below has take from:
% https://stackoverflow.com/questions/20131904/check-if-all-numbers-in-a-list-are-different-in-prolog
all_diff(L) :- \+ (append(_,[X|R],L), memberchk(X,R)).
% the code below has take from:
%   https://stackoverflow.com/questions/4280986/how-to-transpose-a-matrix-in-prolog
transpose([[]|_],[]).
transpose(Matrix, [Row|Rows]):-
    transpose_1st_col(Matrix, Row, RestMatrix),
    transpose(RestMatrix, Rows).
transpose_1st_col([], [], []).
transpose_1st_col([[H|T]|Rows], [H|Hs], [T|Ts]) :- transpose_1st_col(Rows, Hs, Ts).

initiate_board([],_N).
initiate_board([Row|Rows],N):-
    length(Row,N),
    initiate_board(Rows,N).

blocks_to_rows(Board,N,SqrtN,Rows):-   
    length(Rows,N),
    initiate_board(Rows,N),
    loop_a(Rows,0,N,SqrtN,Board).
loop_a(_Rows,N,N,_SqrtN,_Board).
loop_a(Rows,Count,N,SqrtN,Board):-
    Count<N,
    N_count is Count+SqrtN,
    loop_b(Rows,0,N,SqrtN,Board,Count),
    loop_a(Rows,N_count,N,SqrtN,Board).
loop_b(_Rows,N,N,_SqrtN,_Board,_Count_a).
loop_b(Rows,Count,N,SqrtN,Board,Count_a):-
    Count<N,
    N_count is Count+SqrtN,
    loop_c(Rows,0,SqrtN,Board,Count_a,Count),
    loop_b(Rows,N_count,N,SqrtN,Board,Count_a).
loop_c(_Rows,SqrtN,SqrtN,_Board,_Count_a,_Count_b).
loop_c(Rows,Count,SqrtN,Board,Count_a,Count_b):-
    Count<SqrtN,
    N_count is Count+1,
    loop_d(Rows,0,SqrtN,Board,Count_a,Count_b,Count),
    loop_c(Rows,N_count,SqrtN,Board,Count_a,Count_b).
loop_d(_Rows,SqrtN,SqrtN,_Board,_Count_a,_Count_b,_Count_c).
loop_d(Rows,Count,SqrtN,Board,Count_a,Count_b,Count_c):-
    Count<SqrtN,
    N_count is Count+1,
    Rows_i is Count_b//SqrtN+Count_a,
    Rows_j is Count_c*SqrtN+Count,
    Board_i is Count_a+Count_c,
    Board_j is Count_b+Count,
    nth0(Board_i,Board,Row_board),
    nth0(Board_j,Row_board,S),
    (nonvar(S)->
    nth0(Rows_i,Rows,Row_rows),
    nth0(Rows_j,Row_rows,S),
    loop_d(Rows,N_count,SqrtN,Board,Count_a,Count_b,Count);
    loop_d(Rows,N_count,SqrtN,Board,Count_a,Count_b,Count)).

% Task3 Verify Killer 
verify_killer(killer(Instance),Solution,Verified):-
    subset(Instance,Solution),
    length(Board,9),
    initiate_board(Board,9),
    put_hints(Solution,Board),
    once(is_sol(Board,Verified)),
    (var(Verified)-> Verified = killer; true).

put_hints([],_Board).
put_hints([Hint|Hints],Board):-
    Hint = (cell(I,J)=S),
    nth1(I,Board,Row),
    nth1(J,Row,(I,J,S)),
    put_hints(Hints,Board).

is_sol(Board,Verified):-
    distincts_rows(Board,Verified),
    transpose(Board,Cols),
    distincts_rows(Cols,Verified),
    blocks_to_rows(Board,9,3,Boxes),
    distincts_rows(Boxes,Verified),
    check_moves_loop1(Board,1,9,Verified).

distincts_rows([],_Verified).
distincts_rows([Row|Rows],Verified):-
    (var(Verified)->
    (all_diff(Row)->
    distincts_rows(Rows,Verified);
    Row=[A|As],
    find_dup(A,As,Verified));true).
    
% find_dup(_,[],_):- false.
find_dup(First,List,Verified):- %TODO
    find_dup_loop(First,List,Verified),
    (var(Verified)-> Rest=[Second|Rest], find_dup(Second,Rest,Verified)).

find_dup_loop(_First,[],_Verified).
find_dup_loop(First,[Second|Rest],Verified):-
    (var(Verified)->
    First = (I1,J1,Val1),
    Second = (I2,J2,Val2),
    (Val1=Val2->
    Verified = [cell(I1,J1)=Val1,cell(I2,J2)=Val2];
    find_dup_loop(First,Rest,Verified));true).

check_moves_loop1(_Board,Size,Size,_Verified).
check_moves_loop1(Board,I,Size,Verified):-
    (var(Verified)->
    I<Size,
    N_i is I+1,
    check_moves_loop2(Board,I,1,Size,Verified),
    check_moves_loop1(Board,N_i,Size,Verified);
    true).

check_moves_loop2(_Board,_I,Size,Size,_Verified).
check_moves_loop2(Board,I,J,Size,Verified):-
    (var(Verified)->
    J<Size,
    N_j is J+1,
    nth1(I,Board,Row),
    nth1(J,Row,(I,J,S)),
    check_knight_moves(Board,I,J,S,Verified),
    check_king_moves(Board,I,J,S,Verified),
    check_diff(Board,I,J,S,Verified),
    check_moves_loop2(Board,I,N_j,Size,Verified);
    true).

check_knight_moves(Board,I,J,S,Verified):-
    %Down_right_loc
    (var(Verified)->
    A1 is I+2,
    B1 is J+1,
    (A1<1;A1>9;B1<1;B1>9->
    true;
    nth1(A1,Board,Row1),
    nth1(B1,Row1,DR),
    (DR\=S-> true; Verified=[cell(I,J)=S, cell(A1,B1)=DR]));true),
    (var(Verified)->
    % Down_left_loc
    A2 is I+2,
    B2 is J-1,
    (A2<1;A2>9;B2<1;B2>9->
    true;
    nth1(A2,Board,Row2),
    nth1(B2,Row2,DL),
    (DL\=S-> true; Verified=[cell(I,J)=S, cell(A2,B2)=DL]));true),

    % Right_down_loc
    (var(Verified)->
    A6 is I+1,
    B6 is J+2,
    (A6<1;A6>9;B6<1;B6>9->
    true;
    nth1(A6,Board,Row6),
    nth1(B6,Row6,RD),
    (RD\=S-> true; Verified=[cell(I,J)=S, cell(A6,B6)=RD]));true),
    
    % Left_down_loc
    (var(Verified)->
    A8 is I+1,
    B8 is J-2,
    (A8<1;A8>9;B8<1;B8>9->
    true;
    nth1(A8,Board,Row8),
    nth1(B8,Row8,LD),
    (LD\=S-> true; Verified=[cell(I,J)=S, cell(A8,B8)=LD]));true).

        % Up_right_loc
    % (var(Verified)->
    % A3 is I-2,
    % B3 is J+1,
    % (A3<1;A3>9;B3<1;B3>9->
    % true;
    % nth1(A3,Board,Row3),
    % nth1(B3,Row3,UR),
    % (UR\=S-> true; Verified=[cell(I,J)=S, cell(A3,B3)=UR]));true),
    % Up_left_loc
    % (var(Verified)->
    % A4 is I-2,
    % B4 is J-1,
    % (A4<1;A4>9;B4<1;B4>9->
    % true;
    % nth1(A4,Board,Row4),
    % nth1(B4,Row4,UL),
    % (UL\=S-> true; Verified=[cell(I,J)=S, cell(A4,B4)=UL]));true),
    % Right_up_loc
    % (var(Verified)->
    % A5 is I-1,
    % B5 is J+2,
    % (A5<1;A5>9;B5<1;B5>9->
    % true;
    % nth1(A5,Board,Row5),
    % nth1(B5,Row5,RU),
    % (RU\=S-> true; Verified=[cell(I,J)=S, cell(A5,B5)=RU]));true),
    % Left_up_loc
    % (var(Verified)->
    % A7 is I-1,
    % B7 is J-2,
    % (A7<1;A7>9;B7<1;B7>9->
    % true;
    % nth1(A7,Board,Row7),
    % nth1(B7,Row7,LU),
    % (LU\=S-> true; Verified=[cell(I,J)=S, cell(A7,B7)=LU]));true),

check_king_moves(Board,I,J,S,Verified):-
    % Right_down_loc
    (var(Verified)->
    A6 is I+1,
    B6 is J+1,
    (A6<1;A6>9;B6<1;B6>9->
    true;
    nth1(A6,Board,Row6),
    nth1(B6,Row6,RD),
    (RD\=S-> true; Verified=[cell(I,J)=S, cell(A6,B6)=RD]));true),
    % Left_down_loc
    (var(Verified)->
    A8 is I+1,
    B8 is J-1,
    (A8<1;A8>9;B8<1;B8>9->
    true;
    nth1(A8,Board,Row8),
    nth1(B8,Row8,LD),
    (LD\=S-> true; Verified=[cell(I,J)=S, cell(A8,B8)=LD]));true).

    %Down_loc
    % (var(Verified)->
    % A1 is I+1,
    % B1 is J,
    % (A1<1;A1>9;B1<1;B1>9->
    % true;
    % nth1(A1,Board,Row1),
    % nth1(B1,Row1,D),
    % (D\=S-> true; Verified=[cell(I,J)=S, cell(A1,B2)=D]));true),
    % % Left_loc
    % (var(Verified)->
    % A2 is I,
    % B2 is J-1,
    % (A2<1;A2>9;B2<1;B2>9->
    % true;
    % nth1(A2,Board,Row2),
    % nth1(B2,Row2,L),
    % (L\=S-> true; Verified=[cell(I,J)=S, cell(A2,B2)=L]));true),
    % % Right_loc
    % (var(Verified)->
    % A3 is I,
    % B3 is J+1,
    % (A3<1;A3>9;B3<1;B3>9->
    % true;
    % nth1(A3,Board,Row3),
    % nth1(B3,Row3,R),
    % (R\=S-> true; Verified=[cell(I,J)=S, cell(A3,B3)=R]));true),
    % % Up_loc
    % (var(Verified)->
    % A4 is I-1,
    % B4 is J,
    % (A4<1;A4>9;B4<1;B4>9->
    % true;
    % nth1(A4,Board,Row4),
    % nth1(B4,Row4,Up),
    % (Up\=S-> true; Verified=[cell(I,J)=S, cell(A4,B4)=Up]));true),

    % Left_up_loc
    % (var(Verified)->
    % A7 is I-1,
    % B7 is J-1,
    % (A7<1;A7>9;B7<1;B7>9->
    % true;
    % nth1(A7,Board,Row7),
    % nth1(B7,Row7,LU),
    % (LU\=S-> true; Verified=[cell(I,J)=S, cell(A7,B7)=LU]));true),

    % Right_up_loc
    % (var(Verified)->
    % A5 is I-1,
    % B5 is J+1,
    % (A5<1;A5>9;B5<1;B5>9->
    % true;
    % nth1(A5,Board,Row5),
    % nth1(B5,Row5,RU),
    % (RU\=S-> true; Verified=[cell(I,J)=S, cell(A5,B5)=RU]));true),
    
check_diff(Board,I,J,S,Verified):-
    %Down_loc
    (var(Verified)->
    A1 is I+1,
    B1 is J,
    (A1<1;A1>9;B1<1;B1>9->
    true;
    nth1(A1,Board,Row1),
    nth1(B1,Row1,D),
    Diff1 is abs(D-S),
    (Diff1>1->true; Verified=[cell(I,J)=S, cell(A1,B1)=Diff1]));true),
    % Left_loc
    (var(Verified)->
    A2 is I,
    B2 is J-1,
    (A2<1;A2>9;B2<1;B2>9->
    true;
    nth1(A2,Board,Row2),
    nth1(B2,Row2,L),
    Diff2 is abs(L-S),
    (Diff2>1->true; Verified=[cell(I,J)=S, cell(A2,B2)=Diff2]));true),
    % Right_loc
    (var(Verified)->
    A3 is I,
    B3 is J+1,
    (A3<1;A3>9;B3<1;B3>9->
    true;
    nth1(A3,Board,Row3),
    nth1(B3,Row3,R),
    Diff3 is abs(R-S),
    (Diff3>1->true; Verified=[cell(I,J)=S, cell(A3,B3)=Diff3]));true),
    % Up_loc
    (var(Verified)->
    A4 is I-1,
    B4 is J,
    (A4<1;A4>9;B4<1;B4>9->
    true;
    nth1(A4,Board,Row4),
    nth1(B4,Row4,Up),
    Diff4 is abs(Up-S),
    (Diff4>1->true; Verified=[cell(I,J)=S, cell(A4,B4)=Diff4]));true).

map_loop1(Row,[Rowb],I,I):- length(Rowb,9), map_loop2(Row,Rowb,I,1,I).
map_loop1(Map,[Rowb|Restb],I,Length):-
    I<Length,
    N_i is I+1,
    length(Row,9),
    length(Rowb,9),
    append(Row,Rows,Map),
    % append(Row2,Rows2,Board),
    map_loop2(Row,Rowb,I,1,Length),
    map_loop1(Rows,Restb,N_i,Length).
map_loop2([Cell],[Cellb],I,J,J):- Cell= (cell(I,J)= Value), Cellb=Value.
map_loop2([Cell|Row],[Cellb|Rowb],I,J,Length):-
    J<Length,
    N_j is J+1,
    Cell = (cell(I,J)= Value),
    Cellb = Value, 
    map_loop2(Row,Rowb,I,N_j,Length).

put_hints_CNF([],_Map,_Board).
put_hints_CNF([Hint|Hints],Map,Board):-
    Hint = (cell(I,J)=S),
    Loc is (I-1)*9+J,
    convert_form(S,Val),
    nth1(Loc,Map, cell(I,J) = Val),
    nth1(I,Board,Row),
    nth1(J,Row,Val),
    % nth1(J,Row,S),
    put_hints_CNF(Hints,Map,Board).
convert_form(S,Val):-
    length(Val,9),
    nth1(S,Val,1),
    fill_spaces(Val).

fill_spaces([]).
fill_spaces([A|As]):-
    (var(A)-> A= -1; true),
    fill_spaces(As).

%encode_killer(+,-,-)
encode_killer(Instance, Map, CNF):-
    Instance = killer(Hints),
    length(Map,81),
    length(Board,9),
    % initiate_board(Board,9),
    map_loop1(Map,Board,1,9),
    put_hints_CNF(Hints,Map,Board),
    cnf(Board,CNF).

cnf(Board,CNF):-
    cnf_rows_loop(Board,1,CNF1),
    transpose(Board,Cols),
    cnf_rows_loop(Cols,1,CNF2),
    append(CNF1,CNF2,CNF3),
    blocks_to_rows(Board,9,3,Blocks),
    cnf_rows_loop(Blocks,1,CNF4),
    append(CNF3,CNF4,CNF5),
    killer_moves_loop1(Board,1,9,CNF6),
    append(CNF5,CNF6,CNF).

cnf_rows_loop([Row],9,CNF):- take_nth_bit(Row,9,List), exactly_one(List,CNF).
cnf_rows_loop([Row|Rows],Count,CNF):-
    Count<9,
    take_nth_bit(Row,Count,List),
    exactly_one(List,CNF1),
    N_count is Count+1,
    cnf_rows_loop(Rows,N_count,CNF2),
    append(CNF1,CNF2,CNF).

killer_moves_loop1(Board,Size,Size,CNF):-killer_moves_loop2(Board,Size,1,Size,CNF).
killer_moves_loop1(Board,I,Size,CNF):-
    I<Size,
    N_i is I+1,
    killer_moves_loop2(Board,I,1,Size,CNF1),
    killer_moves_loop1(Board,N_i,Size,CNF2),
    append(CNF1,CNF2,CNF).

killer_moves_loop2(_Board,_I,Size,Size,_CNF).
killer_moves_loop2(Board,I,J,Size,CNF):-
    J<Size,
    N_j is J+1,
    nth1(I,Board,Row),
    nth1(J,Row,S),
    knight_moves_cnf(Board,I,J,S,CNF1),
    king_moves_cnf(Board,I,J,S,CNF2),
    check_diff_cnf(Board,I,J,S,CNF4),
    killer_moves_loop2(Board,I,N_j,Size,CNF6),
    append(CNF1,CNF2,CNF3),
    append(CNF3,CNF4,CNF5),
    append(CNF5,CNF6,CNF).

knight_moves_cnf(Board,I,J,Val,CNF):-
    %Down_right_loc
    A1 is I+2,
    B1 is J+1,
    (A1<1;A1>9;B1<1;B1>9->
    CNF1=[];
    nth1(A1,Board,Row1),
    nth1(B1,Row1,DR),
    compare_numbers(Val,DR,CNF1);true),
    % Down_left_loc
    A2 is I+2,
    B2 is J-1,
    (A2<1;A2>9;B2<1;B2>9->
    CNF2=[];
    nth1(A2,Board,Row2),
    nth1(B2,Row2,DL),
    compare_numbers(Val,DL,CNF2);true),
    % Right_down_loc
    A6 is I+1,
    B6 is J+2,
    (A6<1;A6>9;B6<1;B6>9->
    CNF3=[];
    nth1(A6,Board,Row6),
    nth1(B6,Row6,RD),
    compare_numbers(Val,RD,CNF3);true),
    % Left_down_loc
    A8 is I+1,
    B8 is J-2,
    (A8<1;A8>9;B8<1;B8>9->
    CNF4=[];
    nth1(A8,Board,Row8),
    nth1(B8,Row8,LD),
    compare_numbers(Val,LD,CNF4);true),
    %combining formulas
    append(CNF1,CNF2,CNF5),
    append(CNF5,CNF3,CNF6),
    append(CNF6,CNF4,CNF).

king_moves_cnf(Board,I,J,Val,CNF):-
    % Right_down_loc
    A6 is I+1,
    B6 is J+1,
    (A6<1;A6>9;B6<1;B6>9->
    CNF1=[];
    nth1(A6,Board,Row6),
    nth1(B6,Row6,RD),
    compare_numbers(Val,RD,CNF1);true),
    % Left_down_loc
    A8 is I+1,
    B8 is J-1,
    (A8<1;A8>9;B8<1;B8>9->
    CNF2=[];
    nth1(A8,Board,Row8),
    nth1(B8,Row8,LD),
    compare_numbers(Val,LD,CNF2);true),
    append(CNF1,CNF2,CNF).

check_diff_cnf(Board,I,J,Val,CNF):-
    %Down_loc
    A1 is I+1,
    B1 is J,
    (A1<1;A1>9;B1<1;B1>9->
    true;
    nth1(A1,Board,Row1),
    nth1(B1,Row1,D),
    not_consecutives(Val,D,CNF1);true),
    % Right_loc
    A2 is I,
    B2 is J+1,
    (A2<1;A2>9;B2<1;B2>9->
    true;
    nth1(A2,Board,Row3),
    nth1(B2,Row3,R),
    not_consecutives(Val,R,CNF2);true),
    append(CNF1,CNF2,CNF).
    
    % % Left_loc
    % (var(Verified)->
    % A2 is I,
    % B2 is J-1,
    % (A2<1;A2>9;B2<1;B2>9->
    % true;
    % nth1(A2,Board,Row2),
    % nth1(B2,Row2,L),
    % Diff2 is abs(L-S),
    % (;true),
    % Right_loc
    % % Up_loc
    % (var(Verified)->
    % A4 is I-1,
    % B4 is J,
    % (A4<1;A4>9;B4<1;B4>9->
    % true;
    % nth1(A4,Board,Row4),
    % nth1(B4,Row4,Up),
    % Diff4 is abs(Up-S),
    % (Diff4>1->true; Verified=[cell(I,J)=S, cell(A4,B4)=Diff4]));true).

not_consecutives([A1,A2|As],[B1,B2|Bs],CNF):-
    CNF1=[[-A1,-B2],[-A2,-B1]],
    not_consecutives([A2|As],[B2|Bs],CNF2),
    append(CNF1,CNF2,CNF).
exactly_one(List,Formula):-
    % atleast_one(List,CNF1),
    atmost_one(List,CNF2),
    append([List],CNF2,Formula).
% atleast_one(List,[List]).
% atmost_one([B],[]).                 % last cell doesn't matter
atmost_one([A|As],CNF):-
    length(As,N),
    N>0,
    (N>1->
    atmost_one_loop(A,As,CNF1),
    atmost_one(As,CNF2),
    append(CNF1,CNF2,CNF);atmost_one_loop(A,As,CNF)).

atmost_one_loop(A,[B|Bs],CNF):-
    CNF1 = [[-A,-B]], 
    (Bs\=[]->
    atmost_one_loop(A,Bs,CNF2),
    append(CNF1,CNF2,CNF);CNF=CNF1).

compare_numbers([Bit1],[Bit2],[[-Bit1,-Bit2]]).
compare_numbers([Bit1|Bits1],[Bit2|Bits2],CNF):-
    CNF1 =[[-Bit1,-Bit2]],
    compare_numbers(Bits1,Bits2,CNF2),
    append(CNF1,CNF2,CNF).

take_nth_bit([Row],N,List):- nth1(N,Row,Element), List=[Element].
take_nth_bit([Row|Rows],N,List):-
    nth1(N,Row,Element),
    take_nth_bit(Rows,N,List1),
    append([Element],List1,List).


