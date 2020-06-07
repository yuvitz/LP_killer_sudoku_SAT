user:file_search_path(sat, '/home/yuval/Documents/PL/satsolver').
:- use_module(sat(satsolver)).
% bunch of helper predicates

% the code below has take from:
% https://stackoverflow.com/questions/39435709/how-to-remove-duplicates-from-a-list-in-swi-prolog
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
% the code below has take from:
% https://stackoverflow.com/questions/20131904/check-if-all-numbers-in-a-list-are-different-in-prolog
all_diff(L) :- \+ (append(_,[X|R],L), memberchk(X,R)).
% the code below has take from:
% https://stackoverflow.com/questions/4280986/how-to-transpose-a-matrix-in-prolog
transpose([[]|_],[]).
transpose(Matrix, [Row|Rows]):-
    transpose_1st_col(Matrix, Row, RestMatrix),
    transpose(RestMatrix, Rows).
transpose_1st_col([], [], []).
transpose_1st_col([[H|T]|Rows], [H|Hs], [T|Ts]) :- transpose_1st_col(Rows, Hs, Ts).
% the code below has take from:
% https://stackoverflow.com/questions/30464504/how-to-find-the-nth-element-of-a-list-in-prolog
match([H|_],0,H) :-
    !.
match([_|T],N,H) :-
    N > 0, %add for loop prevention
    N1 is N-1,
    match(T,N1,H).
% %% Task1 %%
% %sudoku_propagate(+,-)
sudoku_propagate(Instance, List):-
    sudoku_propagate_help(Instance,[],Lists,[cell(1,1)=_],0,_ExpC),
    append(Lists,List).

sudoku_propagate_help(Instance,PrevL,Lists,[cell(0,0)=A],Found,ExpC):-
    (Found>0->
    sudoku_propagate_help(Instance,PrevL,Lists,[cell(1,1)=A],0,ExpC);
    Found=:=0->Lists=[],ExpC=[]).
 
sudoku_propagate_help(sudoku(SqrtN,Instance),PrevL,Lists,Loc,Found,ExpC):-
    (Loc\= [cell(0,0)=_],
    append(Instance,PrevL,AllList),
    \+is_checked(AllList,Loc),
    constraints(AllList,SqrtN,Loc,Constraints_list,ExplainCons),!,
    remove_duplicates(Constraints_list,Constraints),!,
    length(Constraints,N_curr),
    N is SqrtN**2,
    N_1 is N-1,
    find_next_loc(AllList,Loc,NextLoc,N),
    N_1==N_curr)->                                              
    (find_assignment(Constraints,N,Loc),
    ExpC=[ExplainCons->Loc|ResExp],
    NFound is Found+1,
    append(PrevL,Loc,NewPrev),
    Lists=[Loc|ResList],
    sudoku_propagate_help(sudoku(SqrtN,Instance),NewPrev,ResList,NextLoc,NFound,ResExp));
    (Loc\= [cell(0,0)=_],
    append(Instance,PrevL,AllList),
    N is SqrtN**2,
    find_next_loc(AllList,Loc,NextLoc,N),
    sudoku_propagate_help(sudoku(SqrtN,Instance),PrevL,Lists,NextLoc,Found,ExpC)).

%find next free space
find_next_loc(AllList,[cell(I1,J1)=A1],Loc,N):-
    J1=:=N,I1=:=N->Loc=[cell(0,0)=_B];         
    J1=:=N-> NextJ is 1,
    NextI is I1+1,
    (\+is_checked(AllList,cell(NextI,NextJ)=_)-> Loc=[cell(NextI,NextJ)=_];
    (find_next_loc(AllList,[cell(NextI,NextJ)=A1],Loc,N)));   
    J1<N-> Next is J1+1,
    (\+is_checked(AllList,[cell(I1,Next)=A1])-> Loc=[cell(I1,Next)=_C];
    (find_next_loc(AllList,[cell(I1,Next)=_],Loc,N))).
    
find_assignment(Constraints,N,[cell(_I2,_J2)=A2]):-
    numlist(1,N,L),
    delete_elements(Constraints,L,A2).

delete_elements([X|Constraints],L,A2):-
    delete(L,X,Res),
    delete_elements(Constraints,Res,A2).
delete_elements([],[A2],A2).    

%find the constraints for the current cell   
constraints([cell(I1,J1)=A1|Res1],SqrtN,[cell(I2,J2)=A2],[A1|Constraints_list],[cell(I1,J1)=A1|ExplainCons]):-
    (I1==I2;J1==J2;find_box(I1,J1,I2,J2,SqrtN)),           
    constraints(Res1,SqrtN,[cell(I2,J2)=A2],Constraints_list,ExplainCons).

constraints([],_SqrtN,_Loc,[],[]).
constraints([cell(I1,J1)=_A1|Res1],SqrtN,[cell(I2,J2)=A2],Constraints_list,ExplainCons):-
    (I1\=I2,J1\=J2,\+find_box(I1,J1,I2,J2,SqrtN)),        
    constraints(Res1,SqrtN,[cell(I2,J2)=A2],Constraints_list,ExplainCons).

find_box(I1,J1,RowP,ColP,SqrtN):-
    find_box_helper(1,SqrtN,RowP,SqrtN,Row_box),
    find_box_helper(1,SqrtN,ColP,SqrtN,Col_box),
    member(I1,Row_box),member(J1,Col_box).

find_box_helper(Start,End,I1,SqrtN,Row_box):-
    (numlist(Start,End,L1),member(I1,L1))->Row_box=L1;
    NewStart is Start+SqrtN,
    NewEnd is End+SqrtN,
    find_box_helper(NewStart,NewEnd,I1,SqrtN,Row_box). 

is_checked([],[_Loc]):-false.
is_checked([cell(I1,J1)=_A1|_Res1],[cell(I2,J2)=_A2]):-
    I1==I2,J1==J2.
is_checked([cell(I1,J1)=_A1|Res1],[cell(I2,J2)=A2]):-
        (I1\=I2;J1\=J2) -> is_checked(Res1,[cell(I2,J2)=A2]).
    
%*****************sudoku_propagate_explain(Instance, Explain)***********************
sudoku_propagate_explain(sudoku(SqrtN, Hints), Explain):-
    sudoku_propagate_help(sudoku(SqrtN, Hints),[],_List1,[cell(1,1)=_A],0,Exp),
    N is SqrtN**2,
    no_dups_explain(Exp,Explain,N).

%remove duplcates from the explain
no_dups_explain([],[],_N).
no_dups_explain([Constraints->Loc|Res1],No_dup_exp,N):-
    length(Constraints,Len),
    Len\=N-1,numlist(1,N,NumList),no_dups(Constraints,New_ConsList,NumList),
    No_dup_exp=[New_ConsList->Loc|Res2],no_dups_explain(Res1,Res2,N). 

no_dups_explain([Constraints->Loc|Res1],No_dup_exp,N):-
    length(Constraints,Len),
    Len==N-1,No_dup_exp=[Constraints->Loc|Res2],no_dups_explain(Res1,Res2,N).

no_dups([],[],_).
no_dups([cell(P1,P2)=A|Constraints],New_ConsList,NumList):-
    (member(A,NumList)->
    delete(NumList,A,NewNumlist),
    New_ConsList =[cell(P1,P2)=A|Res],
    no_dups(Constraints,Res,NewNumlist));
    \+member(A,NumList)->no_dups(Constraints,New_ConsList,NumList).    

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
    loop_d(Rows,N_count,SqrtN,Board,Count_a,Count_b,Count_c);
    loop_d(Rows,N_count,SqrtN,Board,Count_a,Count_b,Count_c)).

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
    (take_vals(Row,Row_vals),
    all_diff(Row_vals)->
    distincts_rows(Rows,Verified);
    Row=[A|As],
    find_dup(A,As,Verified));true).

take_vals([A|As],Row_vals):-
    A=(_,_,Val),
    (As\=[]->
    take_vals(As,Rest_vals),
    append([Val],Rest_vals,Row_vals);
    Row_vals=[Val]).
% find_dup(_,[],_):- false.
find_dup(First,List,Verified):- 
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
    nth1(B1,Row1,(_,_,DR)),
    (DR\=S-> true; Verified=[cell(I,J)=S, cell(A1,B1)=DR]));true),
    (var(Verified)->
    % Down_left_loc
    A2 is I+2,
    B2 is J-1,
    (A2<1;A2>9;B2<1;B2>9->
    true;
    nth1(A2,Board,Row2),
    nth1(B2,Row2,(_,_,DL)),
    (DL\=S-> true; Verified=[cell(I,J)=S, cell(A2,B2)=DL]));true),
    % Right_down_loc
    (var(Verified)->
    A6 is I+1,
    B6 is J+2,
    (A6<1;A6>9;B6<1;B6>9->
    true;
    nth1(A6,Board,Row6),
    nth1(B6,Row6,(_,_,RD)),
    (RD\=S-> true; Verified=[cell(I,J)=S, cell(A6,B6)=RD]));true),
    % Left_down_loc
    (var(Verified)->
    A8 is I+1,
    B8 is J-2,
    (A8<1;A8>9;B8<1;B8>9->
    true;
    nth1(A8,Board,Row8),
    nth1(B8,Row8,(_,_,LD)),
    (LD\=S-> true; Verified=[cell(I,J)=S, cell(A8,B8)=LD]));true).

check_king_moves(Board,I,J,S,Verified):-
    % Right_down_loc
    (var(Verified)->
    A6 is I+1,
    B6 is J+1,
    (A6<1;A6>9;B6<1;B6>9->
    true;
    nth1(A6,Board,Row6),
    nth1(B6,Row6,(_,_,RD)),
    (RD\=S-> true; Verified=[cell(I,J)=S, cell(A6,B6)=RD]));true),
    % Left_down_loc
    (var(Verified)->
    A8 is I+1,
    B8 is J-1,
    (A8<1;A8>9;B8<1;B8>9->
    true;
    nth1(A8,Board,Row8),
    nth1(B8,Row8,(_,_,LD)),
    (LD\=S-> true; Verified=[cell(I,J)=S, cell(A8,B8)=LD]));true).

check_diff(Board,I,J,S,Verified):-
    %Down_loc
    (var(Verified)->
    A1 is I+1,
    B1 is J,
    (A1<1;A1>9;B1<1;B1>9->
    true;
    nth1(A1,Board,Row1),
    nth1(B1,Row1,(_,_,D)),
    Diff1 is abs(D-S),
    (Diff1>1->true; Verified=[cell(I,J)=S, cell(A1,B1)=Diff1]));true),
    % Left_loc
    (var(Verified)->
    A2 is I,
    B2 is J-1,
    (A2<1;A2>9;B2<1;B2>9->
    true;
    nth1(A2,Board,Row2),
    nth1(B2,Row2,(_,_,L)),
    Diff2 is abs(L-S),
    (Diff2>1->true; Verified=[cell(I,J)=S, cell(A2,B2)=Diff2]));true),
    % Right_loc
    (var(Verified)->
    A3 is I,
    B3 is J+1,
    (A3<1;A3>9;B3<1;B3>9->
    true;
    nth1(A3,Board,Row3),
    nth1(B3,Row3,(_,_,R)),
    Diff3 is abs(R-S),
    (Diff3>1->true; Verified=[cell(I,J)=S, cell(A3,B3)=Diff3]));true),
    % Up_loc
    (var(Verified)->
    A4 is I-1,
    B4 is J,
    (A4<1;A4>9;B4<1;B4>9->
    true;
    nth1(A4,Board,Row4),
    nth1(B4,Row4,(_,_,Up)),
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
map_loop2([Cell],[Cellb],I,J,J):- 
    length(Value,9), Cell= (cell(I,J)= Value), Cellb=Value.
map_loop2([Cell|Row],[Cellb|Rowb],I,J,Length):-
    J<Length,
    N_j is J+1,
    length(Value,9),
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
    append(CNF5,CNF6,CNF7),
    exactly_one_in_cell(Board,CNF8),
    append(CNF7,CNF8,CNF).

exactly_one_in_cell([],[]).
exactly_one_in_cell([Row|Rows],CNF):-
    exactly_one_in_cell2(Row,CNF1),
    exactly_one_in_cell(Rows,CNF2),
    append(CNF1,CNF2,CNF).
exactly_one_in_cell2([],[]).
exactly_one_in_cell2([Cell|Cells],CNF):-
    exactly_one(Cell,CNF1),
    exactly_one_in_cell2(Cells,CNF2),
    append(CNF1,CNF2,CNF).

cnf_rows_loop([Row],9,CNF):- take_nth_bit(Row,9,List), exactly_one(List,CNF).
cnf_rows_loop([Row],Count,CNF):- 
    Count<9,
    take_nth_bit(Row,Count,List),
    exactly_one(List,CNF1),
    N_count is Count+1,
    cnf_rows_loop([Row],N_count,CNF2),
    append(CNF1,CNF2,CNF).
cnf_rows_loop([Row|Rows],Count,CNF):-
    Count<9->
    take_nth_bit(Row,Count,List),
    exactly_one(List,CNF1),
    N_count is Count+1,
    cnf_rows_loop([Row|Rows],N_count,CNF2),
    append(CNF1,CNF2,CNF);
    take_nth_bit(Row,Count,List),
    exactly_one(List,CNF1),
    cnf_rows_loop(Rows,1,CNF2),
    append(CNF1,CNF2,CNF).

killer_moves_loop1(Board,Size,Size,CNF):-killer_moves_loop2(Board,Size,1,Size,CNF).
killer_moves_loop1(Board,I,Size,CNF):-
    I<Size,
    N_i is I+1,
    killer_moves_loop2(Board,I,1,Size,CNF1),
    killer_moves_loop1(Board,N_i,Size,CNF2),
    append(CNF1,CNF2,CNF).

killer_moves_loop2(_Board,_I,Size,Size,[]).
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
    compare_numbers(Val,DR,CNF1)),
    % Down_left_loc
    A2 is I+2,
    B2 is J-1,
    (A2<1;A2>9;B2<1;B2>9->
    CNF2=[];
    nth1(A2,Board,Row2),
    nth1(B2,Row2,DL),
    compare_numbers(Val,DL,CNF2)),
    % Right_down_loc
    A6 is I+1,
    B6 is J+2,
    (A6<1;A6>9;B6<1;B6>9->
    CNF3=[];
    nth1(A6,Board,Row6),
    nth1(B6,Row6,RD),
    compare_numbers(Val,RD,CNF3)),
    % Left_down_loc
    A8 is I+1,
    B8 is J-2,
    (A8<1;A8>9;B8<1;B8>9->
    CNF4=[];
    nth1(A8,Board,Row8),
    nth1(B8,Row8,LD),
    compare_numbers(Val,LD,CNF4)),
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
    compare_numbers(Val,RD,CNF1)),
    % Left_down_loc
    A8 is I+1,
    B8 is J-1,
    (A8<1;A8>9;B8<1;B8>9->
    CNF2=[];
    nth1(A8,Board,Row8),
    nth1(B8,Row8,LD),
    compare_numbers(Val,LD,CNF2)),
    append(CNF1,CNF2,CNF).

check_diff_cnf(Board,I,J,Val,CNF):-
    %Down_loc
    A1 is I+1,
    B1 is J,
    (A1<1;A1>9;B1<1;B1>9->
    true;
    nth1(A1,Board,Row1),
    nth1(B1,Row1,D),
    not_consecutives(Val,D,CNF1)),
    % Right_loc
    A2 is I,
    B2 is J+1,
    (A2<1;A2>9;B2<1;B2>9->
    true;
    nth1(A2,Board,Row3),
    nth1(B2,Row3,R),
    not_consecutives(Val,R,CNF2)),
    append(CNF1,CNF2,CNF).

not_consecutives([_A],[_B],[]).
not_consecutives([A1,A2|As],[B1,B2|Bs],CNF):-
    CNF1=[[-A1,-B2],[-A2,-B1]],
    not_consecutives([A2|As],[B2|Bs],CNF2),
    append(CNF1,CNF2,CNF).

exactly_one(List,CNF):-
    atmost_one(List,CNF1),
    append([List],CNF1,CNF).
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

%*************************solve_killer(Instance, Solution)************

solve_killer(Instance,Solution) :-
    encode_killer(Instance,Map,CNF),
    sat(CNF),
    decode_killer(Map, Solution),
    verify_killer(Instance, Solution, Verified),Verified = killer.

decode_killer([],[]).
decode_killer([cell(I,J)=A|Res_map], [cell(I,J)=B|Solution]):-
    nth1(B,A,1),
    % ((A==[1,-1,-1,-1,-1,-1,-1,-1,-1])->B=1;
    % (A==[-1,1,-1,-1,-1,-1,-1,-1,-1])->B=2;
    % (A==[-1,-1,1,-1,-1,-1,-1,-1,-1])->B=3;
    % (A==[-1,-1,-1,1,-1,-1,-1,-1,-1])->B=4;
    % (A==[-1,-1,-1,-1,1,-1,-1,-1,-1])->B=5;
    % (A==[-1,-1,-1,-1,-1,1,-1,-1,-1])->B=6;
    % (A==[-1,-1,-1,-1,-1,-1,1,-1,-1])->B=7;
    % (A==[-1,-1,-1,-1,-1,-1,-1,1,-1])->B=8;
    % (A==[-1,-1,-1,-1,-1,-1,-1,-1,1])->B=9),
    decode_killer(Res_map,Solution).

%********************legal_killer(Instance, IsLegal)*****************
legal_killer(Instance, IsLegal):-
    solve_killer(Instance,Sol1),
    encode_killer(Instance,MapSec,CNF),
    add_constraints(MapSec,Sol1,NewCNF),
    append(CNF,[NewCNF],CNF2),
    check(Instance,MapSec,CNF2,Sol1,IsLegal).

%check if the second result is kiiler sudoku constraints
check(Instance,MapSec,CNF2,Sol1,IsLegal):-
        % sat(CNF2),
        (\+solve_killer(Instance,SecSol)->IsLegal=legal,true;
        decode_killer(MapSec,SecSol),
        find_other_assignment(Sol1,SecSol,IsLegal)).

check(_,_MapSec,CNF2,_FirstSolution,IsLegal):-
        % \+(sat(CNF2)),
        IsLegal=legal.    

 %find double assighnment
find_other_assignment([cell(I,J)=Val|Sol1],[cell(I,J)=Sol|SecondSolution],IsLegal):-
    ((Val==Sol)->find_other_assignment(Sol1,SecondSolution,IsLegal),true;
    IsLegal=[cell(I,J)=Val,cell(I,J)=Sol]).

%add constraints that the the solution be diffrent from the first solution
add_constraints([],[],[]).
add_constraints([cell(I,J)=Val|Map],[cell(I,J)=Sol|Sol1],NewCNF):-
    Place_Sol is Sol-1,
    % match(Val,Place_Sol,Ind),
    nth1(Ind,Place_Sol,Val),
    (var(Ind)->NewCNF=[-Ind|ResCnf],add_constraints(Map,Sol1,ResCnf),true;
    add_constraints(Map,Sol1,NewCNF)).

%*********************generate_killer(K, Hints)*********************
generate_killer(K, Hints):-
    gen_hints(K,[],Hints),legal_killer(killer(Hints), IsLegal),IsLegal = legal.

gen_hints(0,_,[]).
gen_hints(K,Prev,[cell(I,J)=Val|Hints]):-
    K>0,
    between(1,9,I),between(1,9,J),between(1,9,Val),
    \+memberchk([I,J],Prev),
    NewK is K-1,
    append(Prev,[[I,J]],Nprev),
    gen_hints(NewK,Nprev,Hints).

lists_firsts_rests([], [], []).
lists_firsts_rests([[F|Os]|Rest], [F|Fs], [Os|Oss]) :-
        lists_firsts_rests(Rest, Fs, Oss).
