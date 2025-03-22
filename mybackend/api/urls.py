from django.urls import path
from rest_framework_simplejwt.views import TokenObtainPairView, TokenRefreshView
from .views import RegisterView, get_smiles_list, get_structure, SmilesQueryView, get_query_history

urlpatterns = [
    path('register/', RegisterView.as_view(), name='register'),
    path('login/', TokenObtainPairView.as_view(), name='login'),
    path('token/refresh/', TokenRefreshView.as_view(), name='token_refresh'),

    # 获取 SMILES 列表（不需要登录）
    path("smiles/", get_smiles_list, name="get_smiles_list"),

    # 解析 SMILES 并返回结构图（不需要登录）
    path("structure/", get_structure, name="get_structure"),

    # 用户输入 SMILES 查询并保存查询记录（需要登录）
    path("query_smiles/", SmilesQueryView.as_view(), name="query_smiles"),

    # 获取当前用户的查询历史（需要登录）
    path("query_history/", get_query_history, name="query_history"),
]
