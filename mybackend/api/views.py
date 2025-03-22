import os
import io
import pandas as pd
from django.conf import settings
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse, HttpResponse
from django.shortcuts import render
from django.utils.decorators import method_decorator
from rdkit import Chem
from rdkit.Chem import Draw
from rest_framework.decorators import api_view, permission_classes
from rest_framework_simplejwt.views import TokenObtainPairView, TokenRefreshView
from django.contrib.auth.models import User
from rest_framework import generics
from rest_framework.permissions import AllowAny, IsAuthenticated
from rest_framework.serializers import ModelSerializer
from rest_framework import status
from rest_framework.response import Response
from rest_framework.views import APIView

from api.models import CustomUser, QueryHistory
from api.serializers import UserSerializer


class UserSerializer(ModelSerializer):
    class Meta:
        model = CustomUser
        fields = ('id', 'username', 'password')
        extra_kwargs = {'password': {'write_only': True}}

class RegisterView(APIView):
    permission_classes = []  # 允许未认证用户访问
    
    def post(self, request):
        serializer = UserSerializer(data=request.data)
        if serializer.is_valid():
            # 创建新用户
            user = CustomUser.objects.create_user(
                username=serializer.validated_data['username'],
                password=serializer.validated_data['password']
            )
            return Response({
                'message': '注册成功',
                'user_id': user.id
            }, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

EXCEL_PATH = os.path.join(settings.BASE_DIR, "data_Chi.csv")

def load_excel_data():
    """读取 Excel 文件并提取去重的 SMILES 前后缀部分"""
    try:
        df = pd.read_csv(EXCEL_PATH)  # 读取整个 CSV 文件
        if "ps_pair" not in df.columns:  # 确保 "ps_pair" 这一列存在
            return set(), set()

        unique_prefixes = set()
        unique_suffixes = set()

        for smiles_string in df["ps_pair"].dropna():
            parts = smiles_string.split("_")
            if len(parts) == 2:
                prefix, suffix = parts

                if prefix:
                    unique_prefixes.add(prefix)
                if suffix:
                    unique_suffixes.add(suffix)

        return unique_prefixes, unique_suffixes
    except Exception as e:
        print(f"读取 Excel 失败: {e}")
        return set(), set()


def check_smiles(smiles):
    """检查 SMILES 字符串的有效性，并返回 RDKit 分子对象"""
    if not isinstance(smiles, str) or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        return mol
    except Exception:
        return None

# API: 获取唯一的前缀和后缀 SMILES 列表
@api_view(["GET"])
def get_smiles_list(request):
    """返回 Excel 表中去重的前缀和后缀 SMILES 数据"""
    unique_prefixes, unique_suffixes = load_excel_data()
    return JsonResponse({"prefixes": list(unique_prefixes), "suffixes": list(unique_suffixes)})

# API: 根据 SMILES 生成化学结构图
@api_view(["GET"])
def get_structure(request):
    """解析 SMILES 并返回化学结构图（不保存查询历史）"""
    smiles_string = request.GET.get("smiles", "").strip()

    if not smiles_string:
        return JsonResponse({"error": "未提供 SMILES"}, status=400)

    mol = check_smiles(smiles_string)
    if mol is None:
        return JsonResponse({"error": "无效的 SMILES"}, status=400)

    # 生成图片
    img_io = io.BytesIO()
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save(img_io, format="PNG")
    img_io.seek(0)

    return HttpResponse(img_io.getvalue(), content_type="image/png")


class SmilesQueryView(APIView):
    permission_classes = [IsAuthenticated]  # 使用 DRF 的认证机制
    
    def get(self, request):
        smiles_string = request.GET.get("smiles", "").strip()

        if not smiles_string:
            return JsonResponse({"error": "未提供 SMILES"}, status=400)

        mol = check_smiles(smiles_string)
        if mol is None:
            return JsonResponse({"error": "无效的 SMILES"}, status=400)

        # 记录查询历史
        QueryHistory.objects.create(user=request.user, smiles=smiles_string)

        # 生成图片
        img_io = io.BytesIO()
        img = Draw.MolToImage(mol, size=(300, 300))
        img.save(img_io, format="PNG")
        img_io.seek(0)

        return HttpResponse(img_io.getvalue(), content_type="image/png")

@api_view(["GET"])
@permission_classes([IsAuthenticated])  # 需要用户登录
def get_query_history(request):
    """获取当前用户的查询历史"""
    history = QueryHistory.objects.filter(user=request.user).order_by("-created_at")

    return Response({
        "history": [
            {"smiles": record.smiles, "created_at": record.created_at.strftime("%Y-%m-%d %H:%M:%S")}
            for record in history
        ]
    })
