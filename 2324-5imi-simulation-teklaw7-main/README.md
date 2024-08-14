[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/3jxxEx1R)
[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-718a45dd9cf7e7f842a935f5ebbe5719a5e09af4491e668f4dbf3b35d5cca122.svg)](https://classroom.github.com/online_ide?assignment_repo_id=12463771&assignment_repo_type=AssignmentRepo)
# 5imi simulation
TP noté de simulation de tissu

## Compilation

Pour compiler et executer à partir du CMakeLists.txt

```sh
mkdir build
cd build
cmake ..
make
cd ..
./build/programme_1
```

ou 

```sh
mkdir build
cmake . -B build
make -C ./build && ./build/programme_1
```

**Note sur l'utilisation des IDE (QtCreator, etc)**

Le repertoire d'execution doit être ici.
C'est a dire que le repertoire data/ doit être accessible.
