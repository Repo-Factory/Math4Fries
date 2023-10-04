snake_to_uppercase() 
{
    local snake_case="$1"
    uppercase_string="${snake_case^^}"
    echo "$uppercase_string"
}

snake_to_camel() 
{
    local snake_case="$1"
    local parts=()
    IFS="_" read -ra parts <<< "$snake_case"
    local camel_case=""
    for part in "${parts[@]}"; do
      camel_case+="$(tr '[:lower:]' '[:upper:]' <<< ${part:0:1})${part:1}"
    done
    echo "$camel_case"
}