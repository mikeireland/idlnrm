pro make_format_string,dummy
;;just change the variable "file" to whatever you want
;; A18,F8.3,
;;RAhour, DEdeg,
;;RAhour:RAhour, 
file = 'rv_format.txt'
;;openr, 1, file
;;nrows = file_lines(file)
;;lines = strarr(nrows)
;;readf, 1, lines
;;close,1
readcol,file,byt,form,unit,vars,desc,format='(A,A,A,A,A)'

format_str = form[0]
var_str = vars[0]
struct_str = vars[0]
structure = vars[0]+':'+form[0]+','


for i = 1,n_elements(form)-1 do begin
    var_str = var_str[0] + ',' + vars[i]
    format_str = format_str + ',' + form[i] 
    struct_str = struct_str[0] + ':' + vars[i]
    structure = structure + vars[i] + ':' +form[i] + ',' 
endfor

print,'formats_str:',format_str
print,'var_str:',var_str
print,'struct_str',struct_str
print,'structure=', structure
end
